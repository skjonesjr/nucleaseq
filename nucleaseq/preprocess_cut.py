import sys
import os
import random
import yaml
import multiprocessing as mp
import numpy as np
import pandas as pd
from collections import Counter, defaultdict
from freebarcodes import editmeasures
from freebarcodes.decode import FreeDivBarcodeDecoder
from nucleaseq.seqtools import dna_rev_comp, bases, load_read_name_seq_items, load_sample_given_read_name
from nucleaseq.OligosContainer import OligosContainer
from nucleaseq.NucleaSeqOligo import NucleaSeqOligo
import logging


log = logging.getLogger()
if not log.handlers:
    handler = logging.StreamHandler()
    formatter = logging.Formatter("%(asctime)s   %(message)s", "%Y-%m-%d %H:%M:%S")
    handler.setFormatter(formatter)
    log.addHandler(handler)
    log.setLevel(logging.INFO)

#-------------------------------------------------------------------------------
# Arguments
usg_fmt = '{} <targets_file> <target_name> <pamtarg_pos> <exploded_oligos_file> <max_primer_err> <max_bc_err> <read_names_by_seq_file> <read_names_by_sample_file> <out_dir_prefix> <start> <inc> <end> <nprocs>'.format(sys.argv[0])
if len(sys.argv) != len(usg_fmt.split()):
    sys.exit('Usage: ' + usg_fmt)

targets_file = sys.argv[1]
target_name = sys.argv[2]
pamtarg_pos = int(sys.argv[3])
exploded_oligos_file = sys.argv[4]
max_primer_err = int(sys.argv[5])
max_bc_err = int(sys.argv[6])
read_names_by_seq_file = sys.argv[7]
read_names_by_sample_file = sys.argv[8]
out_dir_prefix = sys.argv[9]
start, inc, end, nprocs = map(int, sys.argv[10:])
#-------------------------------------------------------------------------------



log.info('NProcs: {}'.format(nprocs or 'All'))

log.info('Loading target details')
targets = yaml.load(open(targets_file))
perfect_target = targets[target_name]
log.info('Target {}: {}'.format(target_name, perfect_target))

def count_alignment_edits(ref_align, observed_align):
    return sum(1 for rc, oc in zip(ref_align, observed_align) if rc != oc)

def identify_side_and_barcode(seq):
    for cr, side in [(oligo_container.cr_left, 'left'),
                     (oligo_container.cr_right, 'right')]:
        res = editmeasures.prefix_identification(cr, seq, max_primer_err)
        if res:
            pos, cr_edits = res
            obs_bc = seq[pos:pos + oligo_container.bc_len]
            bc = bc_decoder.decode(obs_bc.replace('N', 'A'))
            return side, bc, pos, cr_edits
    return None, None, None, None

def left_piece_starts_edits_and_cut_pamtarg_coord(aligned_pieces, oligo):
    # Determine cut piece
    for cut_piece_idx, (_, observed_align) in enumerate(aligned_pieces):
        if set(observed_align) == set('-'):
            cut_piece_idx -= 1
            break
    
    # Separate by type
    full_pieces = aligned_pieces[:cut_piece_idx]
    cut_piece = aligned_pieces[cut_piece_idx]
    gone_pieces = aligned_pieces[cut_piece_idx + 1:]
    
    # Full pieces
    obs_start_idx = 0
    piece_starts, piece_edits = [], []
    for ref_align, observed_align in full_pieces:
        piece_starts.append(obs_start_idx)
        piece_edits.append(count_alignment_edits(ref_align, observed_align))
        obs_start_idx += sum(1 for oc in observed_align if oc != '-')
        
    # Cut piece
    piece_starts.append(obs_start_idx)
    ref_align, observed_align = cut_piece
    cut_idx = 0
    for i, oc in enumerate(observed_align):
        if oc != '-':
            cut_idx = i + 1
    piece_edits.append(count_alignment_edits(ref_align[:cut_idx], observed_align[:cut_idx]))
    
    # Cut pos
    cut_pos = sum(1 for rc in ref_align[:cut_idx] if rc != '-')
    cut_pos += sum(len(p) for p in oligo.pieces[:cut_piece_idx])
    if oligo in oligo_container.target_oligos:
        cut_pamtarg_coord = oligo.convert_oligo_pos_to_pamtarg_coord(cut_pos)
    else:
        cut_pamtarg_coord = cut_pos

    
    # Gone pieces
    for _ in gone_pieces:
        piece_starts.append(np.nan)
        piece_edits.append(np.nan)
    
    assert len(aligned_pieces) == len(piece_starts) == len(piece_edits), ('\n'.join(aligned_pieces), 
                                                                          piece_starts, 
                                                                          piece_edits) 
    return piece_starts[1:], piece_edits, cut_pamtarg_coord  # Don't return that the first piece starts at zero

def right_piece_starts_edits_and_cut_pamtarg_coord(aligned_pieces, oligo):
    # Determine cut piece
    for cut_piece_idx, (_, observed_align) in enumerate(aligned_pieces):
        if set(observed_align) != set('-'):
            break
    
    # Separate by type
    gone_pieces = aligned_pieces[:cut_piece_idx]
    cut_piece = aligned_pieces[cut_piece_idx]
    full_pieces = aligned_pieces[cut_piece_idx + 1:]
    
    # Gone pieces
    piece_starts, piece_edits = [], []
    for _ in gone_pieces:
        piece_starts.append(np.nan)
        piece_edits.append(np.nan)
    
    # Cut piece
    piece_starts.append(0)
    ref_align, observed_align = cut_piece
    cut_idx = 0
    for cut_idx, oc in enumerate(observed_align):
        if oc != '-':
            break
    piece_edits.append(count_alignment_edits(ref_align[cut_idx:], observed_align[cut_idx:]))

    # Cut pos
    cut_pos = sum(1 for rc in ref_align[:cut_idx] if rc != '-')
    cut_pos += sum(len(p) for p in oligo.pieces[:cut_piece_idx])
    if oligo in oligo_container.target_oligos:
        cut_pamtarg_coord = oligo.convert_oligo_pos_to_pamtarg_coord(cut_pos)
    else:
        cut_pamtarg_coord = cut_pos

    # Full pieces
    obs_start_idx = sum(1 for oc in observed_align if oc != '-') # from cut piece
    for ref_align, observed_align in full_pieces:
        piece_starts.append(obs_start_idx)
        piece_edits.append(count_alignment_edits(ref_align, observed_align))
        obs_start_idx += sum(1 for oc in observed_align if oc != '-')
        
    assert len(aligned_pieces) == len(piece_starts) == len(piece_edits), ('\n'.join(aligned_pieces), 
                                                                          piece_starts, 
                                                                          piece_edits) 
    return piece_starts[1:], piece_edits, cut_pamtarg_coord  # Don't return that the first piece starts at zero

def test_oligo_alignments(oligo_container):
    log.debug('Showing test alignments.')
    oligo = oligo_container.oligos[0]
    log.debug('All oligo pieces:\n' + str(oligo.pieces))
    left_oligo = oligo.sequence[:50] + oligo.sequence[55:65]
    right_oligo = oligo.sequence[45:55] + oligo.sequence[60:]

    cut_oligo = left_oligo
    aligned_pieces = oligo.align_and_return_as_pieces(cut_oligo, align_method='global_cfe')
    log.debug('Aligned pieces:\n' + '\n'.join(map(str, aligned_pieces)))
    log.debug('Starts and nEdits:\n' + str(left_piece_starts_edits_and_cut_pamtarg_coord(aligned_pieces, oligo)))

    cut_oligo = right_oligo
    aligned_pieces = oligo.align_and_return_as_pieces(cut_oligo, align_method='global_cfe')
    log.debug('Aligned pieces:\n' + '\n'.join(map(str, aligned_pieces)))
    log.debug('Starts and nEdits:\n' + str(right_piece_starts_edits_and_cut_pamtarg_coord(aligned_pieces, oligo)))


log.info('Loading oligo details')
oligo_container = OligosContainer(exploded_oligos_file,
                                  perfect_target,
                                  pamtarg_pos)
test_oligo_alignments(oligo_container)

log.info('Loading read name properties')
read_name_items = load_read_name_seq_items(read_names_by_seq_file)
sample_given_read_name = load_sample_given_read_name(
    read_names_by_sample_file
)

if start >= len(read_name_items):
    raise ValueError('Start after last seq')

log.info('Loading FREE barcode decoder')
bc_decoder = FreeDivBarcodeDecoder()
bc_decoder.build_codebook_from_codewords(oligo_container.barcodes, max_bc_err)



pieces_names = NucleaSeqOligo.pieces_names

def preprocess_cutseq(tup):
    seq, read_names = tup
    #-----------------------------
    # Process sides
    #-----------------------------
    seq_rc = dna_rev_comp(seq)
    left_side, left_bc, left_pos, left_edits = identify_side_and_barcode(seq)
    right_side, right_bc, right_pos, right_edits = identify_side_and_barcode(seq_rc)
    if set([left_side, right_side]) in [set(['left', 'right']), set(['left']), set(['right'])]:
        if left_side == right_side:
            if left_side == 'left':
                return 'two lefts', len(read_names)
            else:
                return 'two rights', len(read_names)
        return 'full', len(read_names)
        
    if left_side == right_side == None:
        return 'no primers', len(read_names)
        
    #-----------------------------
    # Orient the seq
    #-----------------------------
    if left_side == 'right' or right_side == 'left':
        seq, seq_rc = seq_rc, seq
        left_side, left_bc, left_pos, left_edits, right_side, right_bc, right_pos, right_edits = \
            right_side, right_bc, right_pos, right_edits, left_side, left_bc, left_pos, left_edits

    assert ((left_side == 'left' and right_side is None)
            or (left_side is None and right_side == 'right')), (left_side, right_side, seq)
    
    if left_side == 'left' and left_bc is None:
        return 'left none', len(read_names)
        
    if right_side == 'right' and right_bc is None:
        return 'right none', len(read_names)
        
    #------------------------------------
    # Process and store
    #------------------------------------
    
    read_names = list(read_names)
    samples = [sample_given_read_name[rn] for rn in read_names]

    if left_side == 'left':
        side = 'left'
        oligo = oligo_container.oligo_given_barcode_given_side[left_side][left_bc]
        try:
            aligned_pieces = oligo.align_and_return_as_pieces(seq, align_method='global_cfe')
        except:
            aligned_pieces = oligo.align_and_return_as_pieces(seq, align_method='global')
        piece_starts, piece_edits, cut_pamtarg_coord = left_piece_starts_edits_and_cut_pamtarg_coord(aligned_pieces, 
                                                                                                     oligo)
    else:
        assert right_side == 'right', (left_side, right_side, seq)
        side = 'right'
        oligo = oligo_container.oligo_given_barcode_given_side[right_side][right_bc]
        try:
            aligned_pieces = oligo.align_and_return_as_pieces(seq, align_method='global_cfe')
        except:
            aligned_pieces = oligo.align_and_return_as_pieces(seq, align_method='global')
        piece_starts, piece_edits, cut_pamtarg_coord = right_piece_starts_edits_and_cut_pamtarg_coord(aligned_pieces, 
                                                                                                      oligo)

    edits_names = ['edits_' + pn for pn in pieces_names]
    start_names = ['start_' + pn for pn in pieces_names[1:]]
    aligned_pieces_names = ['aligned_' + pn for pn in pieces_names]
    edits_dict = {name: edits for name, edits in zip(edits_names, piece_edits)}
    starts_dict = {name: start for name, start in zip(start_names, piece_starts)}
    aligned_pieces_dict = {name: [al for _ in range(len(samples))]
                           for name, al in zip(aligned_pieces_names, aligned_pieces)}
        
    columns = [
        'read_name',
        'sample',
        'side',
        'oriented_seq',
        'seq_len',
        'oligo',
        'cut_pamtarg_coord',
        'total_edits',
    ] + edits_names + start_names + aligned_pieces_names
        
    d = {
        'read_name': read_names,
        'sample': samples,
        'side': side,
        'oriented_seq': seq,
        'seq_len': len(seq),
        'oligo': oligo.sequence,
        'cut_pamtarg_coord': cut_pamtarg_coord,
        'total_edits': sum(piece_edits),
    }
    d.update(edits_dict)
    d.update(starts_dict)
    d.update(aligned_pieces_dict)
        
    return side, pd.DataFrame(d, columns=columns)

# Start the races
last_end = end 
log.info('Start idx: {:,d}'.format(start))
log.info('End idx: {:,d}'.format(end))
log.info('Increment: {:,d}'.format(inc))

out_dir = '{}_cut_data_files'.format(out_dir_prefix)
if not os.path.exists(out_dir):
    os.mkdir(out_dir)
seq_idx_digits = len(str(len(read_name_items)))
out_fname_template = 'cut_data.{:0%dd}-{:0%dd}.pkl' % (seq_idx_digits, seq_idx_digits)

while start < last_end:
    end = start + inc
    pl = mp.Pool(nprocs)
    res = pl.map(preprocess_cutseq, read_name_items[start:end])
    pl.close()
    pl.join()
    del pl

    stats = Counter()
    cut_data_frames = []
    for label, df in res:
        if label in ['left', 'right']:
            stats['{} seqs'.format(label)] += 1
            stats['{} reads'.format(label)] += len(df)
            cut_data_frames.append(df)
        elif label in ['full', 'two lefts', 'two rights', 'no primers', 'left none', 'right none']:
            stats['{} seqs'.format(label)] += 1
            stats['{} reads'.format(label)] += df
        else:
            raise ValueError('Unexpected result type "{}"'.format(label))
    stats['cut seqs'] = len(cut_data_frames)
    del res

    if cut_data_frames:
        cut_data = pd.concat(cut_data_frames)
        cut_data.to_pickle(os.path.join(out_dir, out_fname_template.format(start, end-1)))
        del cut_data_frames
        del cut_data
            
    log.info('{:,d}-{:,d}: '.format(start, end) + str(sorted(stats.items())))
    start += inc

