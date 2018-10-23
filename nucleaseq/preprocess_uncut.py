import sys
import os
import random
import yaml
import multiprocessing as mp
import numpy as np
import pandas as pd
from pathos.multiprocessing import ProcessingPool
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
usg_fmt = '{} <start> <inc> <end>'.format(sys.argv[0])
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
start, inc, end, nprocs = map(int, sys.argv[9:])
#-------------------------------------------------------------------------------


log.info('NProcs: {}'.format(nprocs or 'All'))

log.info('Loading target details')
targets = yaml.load(open(targets_file))
perfect_target = targets[target_name]
log.info('Target {}: {}'.format(target_name, perfect_target))

log.info('Loading oligo details')
oligo_container = OligosContainer(exploded_oligos_file,
                                  perfect_target,
                                  pamtarg_pos)

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

def count_alignment_edits(ref_align, observed_align):
    return sum(1 for rc, oc in zip(ref_align, observed_align) if rc != oc)

def piece_starts_and_edits(aligned_pieces):
    obs_start_idx = 0
    piece_starts, piece_edits = [], []
    for ref_align, observed_align in aligned_pieces:
        piece_starts.append(obs_start_idx)
        piece_edits.append(count_alignment_edits(ref_align, observed_align))
        obs_start_idx += sum(1 for oc in observed_align if oc != '-')
    return piece_starts[1:], piece_edits  # Don't return that the first piece starts at zero

pieces_names = NucleaSeqOligo.pieces_names

def process_full_none_wrong(tup):
    seq, read_names = tup
    #-----------------------------
    # Process sides
    #-----------------------------
    seq_rc = dna_rev_comp(seq)
    left_side, left_bc, left_pos, left_edits = identify_side_and_barcode(seq)
    right_side, right_bc, right_pos, right_edits = identify_side_and_barcode(seq_rc)
    if set([left_side, right_side]) != set(['left', 'right']):
        return 'other', len(read_names)
        
    #-----------------------------
    # Orient the seq
    #-----------------------------
    if left_side == 'right':
        seq, seq_rc = seq_rc, seq
        left_side, left_bc, left_pos, left_edits, right_side, right_bc, right_pos, right_edits = \
            right_side, right_bc, right_pos, right_edits, left_side, left_bc, left_pos, left_edits

    #------------------------------------
    # Process and store
    #------------------------------------
    left_oligo = None if left_bc is None else oligo_container.oligo_given_barcode_given_side[left_side][left_bc]
    right_oligo = None if right_bc is None else oligo_container.oligo_given_barcode_given_side[right_side][right_bc]
    
    read_names = list(read_names)
    samples = [sample_given_read_name[rn] for rn in read_names]
    
    if left_oligo and right_oligo and left_oligo == right_oligo:
        # Good sequence
        oligo = left_oligo
        aligned_pieces = oligo.align_and_return_as_pieces(seq)
        piece_starts, piece_edits = piece_starts_and_edits(aligned_pieces)
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
            'oriented_seq',
            'seq_len',
            'oligo',
            'total_edits',
        ] + edits_names + start_names + aligned_pieces_names
        
        d = {
            'read_name': read_names,
            'sample': samples,
            'oriented_seq': seq,
            'seq_len': len(seq),
            'oligo': oligo.sequence,
            'total_edits': sum(piece_edits),
        }
        d.update(edits_dict)
        d.update(starts_dict)
        d.update(aligned_pieces_dict)
        
        return 'full', pd.DataFrame(d, columns=columns)
    if left_oligo is None and right_oligo is None:
        # Both sides bad. Nothing to do.
        return 'both none', len(read_names)
    elif left_oligo is None or right_oligo is None:
        # One side None
        if left_oligo is None:
            oligo = right_oligo
            missing_side = 'left'
        else:
            oligo = left_oligo
            missing_side = 'right'
            
        aligned_pieces = oligo.align_and_return_as_pieces(seq)
        piece_starts, piece_edits = piece_starts_and_edits(aligned_pieces)
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
            'oriented_seq',
            'seq_len',
            'missing_side',
            'oligo',
            'total_edits',
        ] + edits_names + start_names + aligned_pieces_names
        
        d = {
            'read_name': read_names,
            'sample': samples,
            'oriented_seq': seq,
            'seq_len': len(seq),
            'missing_side': missing_side,
            'oligo': oligo.sequence,
            'total_edits': sum(piece_edits),
        }
        d.update(edits_dict)
        d.update(starts_dict)
        d.update(aligned_pieces_dict)
        
        return 'none', pd.DataFrame(d, columns=columns)
    else:
        # One side wrong
        left_aligned_pieces = left_oligo.align_and_return_as_pieces(seq)
        left_piece_starts, left_piece_edits = piece_starts_and_edits(left_aligned_pieces)
        right_aligned_pieces = right_oligo.align_and_return_as_pieces(seq)
        right_piece_starts, right_piece_edits = piece_starts_and_edits(right_aligned_pieces)
        
        left_edits_names = ['left_edits_' + pn for pn in pieces_names]
        left_start_names = ['left_start_' + pn for pn in pieces_names[1:]]
        left_aligned_pieces_names = ['left_aligned_' + pn for pn in pieces_names]
        left_edits_dict = {name: edits for name, edits in zip(left_edits_names, left_piece_edits)}
        left_starts_dict = {name: start for name, start in zip(left_start_names, left_piece_starts)}
        left_aligned_pieces_dict = {name: [al for _ in range(len(samples))]
                                    for name, al in zip(left_aligned_pieces_names, left_aligned_pieces)}
        
        right_edits_names = ['right_edits_' + pn for pn in pieces_names]
        right_start_names = ['right_start_' + pn for pn in pieces_names[1:]]
        right_aligned_pieces_names = ['right_aligned_' + pn for pn in pieces_names]
        right_edits_dict = {name: edits for name, edits in zip(right_edits_names, right_piece_edits)}
        right_starts_dict = {name: start for name, start in zip(right_start_names, right_piece_starts)}
        right_aligned_pieces_dict = {name: [al for _ in range(len(samples))]
                                     for name, al in zip(right_aligned_pieces_names, right_aligned_pieces)}
        
        columns = [
            'read_name',
            'sample',
            'oriented_seq',
            'seq_len',
            'left_oligo',
            'left_total_edits',
        ] + left_edits_names + left_start_names + left_aligned_pieces_names + [
            'right_oligo',
            'right_total_edits'
        ] + right_edits_names + right_start_names + right_aligned_pieces_names
        
        d = {
            'read_name': read_names,
            'sample': samples,
            'oriented_seq': seq,
            'seq_len': len(seq),
            'left_oligo': left_oligo.sequence,
            'left_total_edits': sum(left_piece_edits),
            'right_oligo': right_oligo.sequence,
            'right_total_edits': sum(right_piece_edits),
        }
        d.update(left_edits_dict)
        d.update(left_starts_dict)
        d.update(left_aligned_pieces_dict)
        d.update(right_edits_dict)
        d.update(right_starts_dict)
        d.update(right_aligned_pieces_dict)
        
        return 'wrong', pd.DataFrame(d, columns=columns)

# Start the races
start = start
last_end = end 
log.info('Start idx: {:,d}'.format(start))
log.info('End idx: {:,d}'.format(end))
log.info('Increment: {:,d}'.format(inc))

seq_types = ['full', 'none', 'wrong']
out_dir = {seq_type: '{}_data_files'.format(seq_type) for seq_type in seq_types}
for dname in out_dir.values():
    if not os.path.exists(dname):
        os.mkdir(dname)
seq_idx_digits = len(str(len(read_name_items)))
out_fname_template = {
    seq_type: '%s_data.{:0%dd}-{:0%dd}.pkl' % (seq_type, seq_idx_digits, seq_idx_digits)
    for seq_type in seq_types
}

while start < last_end:
    end = start + inc
    pl = mp.Pool(nprocs)
    #pl = ProcessingPool(nprocs)
    res = pl.map(process_full_none_wrong, read_name_items[start:end])
    pl.close()
    #pl.terminate()
    pl.join()
    del pl

    stats = Counter()
    data_frames = {seq_type: [] for seq_type in seq_types}
    for label, df in res:
        if label in seq_types:
            data_frames[label].append(df)
        elif label in ['both none', 'other']:
            stats['{} seqs'.format(label)] += 1
            stats['{} reads'.format(label)] += df  # df is here an int
        else:
            raise ValueError('Unexpected result type "{}"'.format(label))
    for seq_type in seq_types:
        stats['{} seqs'.format(seq_type)] = len(data_frames[seq_type])

    for seq_type in seq_types:
        if data_frames[seq_type]:
            data = pd.concat(data_frames[seq_type])
            data.to_pickle(os.path.join(out_dir[seq_type],
                                        out_fname_template[seq_type].format(start, end-1)))
            del data
    del data_frames
        
    del res

    log.info('{:,d}-{:,d}: '.format(start, end) + str(sorted(stats.items())))
    start += inc
