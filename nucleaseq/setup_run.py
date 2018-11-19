import sys
import os
import random
import itertools
import pickle
import logging
import numpy as np
from collections import defaultdict
from Bio import SeqIO
import nucleaseq
import misc
from nucleaseq.seqtools import bases
from freebarcodes.editmeasures import simple_hamming_distance


log = logging.getLogger(__name__)

THIS_DIR = os.path.dirname(os.path.abspath(__file__))

def rand_seq(seq_len):
    return ''.join(random.choice(bases) for _ in xrange(seq_len))


def get_max_ham_dists(min_len, max_len):
    dists = defaultdict(list)
    for _ in xrange(50000):
        ref_seq = rand_seq(max_len)
        new_seq = rand_seq(max_len)
        for i in range(min_len, max_len+1):
            dists[i].append(simple_hamming_distance(ref_seq[:i], new_seq[:i]))
    max_ham_dists = [min(np.percentile(dists[i], 0.05), int(i/4)) for i in range(min_len, max_len+1)]
    return max_ham_dists


def isint(a):
    try:
        int(a)
        return float(a) == int(a)
    except:
        return False


def stitch_reads(arguments):
    """
    Classifies reads by overlapping ML sequence identity.
    """
    assert all(map(isint, (arguments.min_len, arguments.max_len))), 'Min and Max lens must be integers'
    min_allowed_overlap = 8
    #--------------------------------------------------------------------------------
    # Load log_p dict of dicts of lists. Addessed as follows:
    #   
    #   log_p_struct[true_base][read_base][phred_score]
    #--------------------------------------------------------------------------------
    log_p_fpath = os.path.join(os.path.abspath(os.path.join(THIS_DIR, '..')), 'resources', 'base_logp.pkl')
    #log_p_fpath = os.path.join(THIS_DIR, 'resources', 'base_logp.pkl')
    with open(log_p_fpath) as f:
        log_p_struct = pickle.load(f)

    #--------------------------------------------------------------------------------
    # Find max hamming distances per length considered
    #--------------------------------------------------------------------------------
    log.info('Finding max hamming distances by length')
    max_ham_dists = get_max_ham_dists(0, arguments.max_len)
    log.info('Max hamming distances by length found')

    #--------------------------------------------------------------------------------
    # Make classifying function with given params
    #--------------------------------------------------------------------------------
    bases = 'ACGT'
    bases_set = set(bases)
    def classify_seq(rec1, rec2):
        # Store as strings
        seq1 = str(rec1.seq)
        seq2_rc = str(rec2.seq.reverse_complement())
        loc_max_len = min(arguments.max_len, len(seq1), len(seq2_rc))

        # Find aligning sequence, indels are not allowed, starts of reads included
        sig_lens = [i for i in xrange(arguments.min_len, loc_max_len + 1)
                    if simple_hamming_distance(seq1[:i], seq2_rc[-i:]) < max_ham_dists[i]]
        if loc_max_len < arguments.max_len:
            # max_len longer than seq1 or seq2
            sig_lens.extend(
                [len(seq1) + len(seq2_rc) - i for i in xrange(min_allowed_overlap, loc_max_len)
                 if simple_hamming_distance(seq1[-i:], seq2_rc[:i]) < max_ham_dists[i]]
            )
        if len(sig_lens) != 1:
            return None

        sig_len = sig_lens[0]
        complete_overlap = bool(sig_len <= loc_max_len) 

        if complete_overlap:
            seq1_match = seq1[:sig_len]
            seq2_match = seq2_rc[-sig_len:]
            quals1 = rec1.letter_annotations['phred_quality'][:sig_len]
            quals2 = rec2.letter_annotations['phred_quality'][::-1][-sig_len:]
        else:
            overlap_len = len(seq1) + len(seq2_rc) - sig_len 
            seq1_match = seq1[-overlap_len:]
            seq2_match = seq2_rc[:overlap_len]
            quals1 = rec1.letter_annotations['phred_quality'][-overlap_len:]
            quals2 = rec2.letter_annotations['phred_quality'][::-1][:overlap_len]

        # Build concensus sequence
        ML_bases = []
        for r1, q1, r2, q2 in zip(seq1_match, quals1, seq2_match, quals2):
            if r1 in bases and r1 == r2:
                ML_bases.append(r1)
            elif set([r1, r2]) <= bases_set and q1 > 2 and q2 > 2:
                r1_score = log_p_struct[r1][r1][q1] + log_p_struct[r1][r2][q2]
                r2_score = log_p_struct[r2][r1][q1] + log_p_struct[r2][r2][q2]
                if r1_score > r2_score:
                    ML_bases.append(r1)
                else:
                    ML_bases.append(r2)
            elif r1 in bases and q1 > 2:
                ML_bases.append(r1)
            elif r2 in bases and q2 > 2:
                ML_bases.append(r2)
            else:
                return None

        if complete_overlap:
            return ''.join(ML_bases)
        else:
            return seq1[:-overlap_len] + ''.join(ML_bases) + seq2_rc[overlap_len:]

    #--------------------------------------------------------------------------------
    # Pair fpaths and classify seqs
    #--------------------------------------------------------------------------------
    fastq_fpaths = [fpath
                    for sample_dir in arguments.sample_dirs
                    for fpath in glob.glob(sample_dir, '*')]

    pe_fpaths, se_fpaths = misc.find_paired_and_unpaired_files_from_fpaths(fastq_fpaths)
    assert not se_fpaths, 'All fastq files should be paired:\n{}'.format('\n'.join(se_fpaths))

    read_names_given_seq = defaultdict(list)
    for fpath1, fpath2 in pe_fpaths:
        log.info('{}, {}'.format(*map(os.path.basename, (fpath1, fpath2))))
        discarded = 0
        total = 0
        for i, (rec1, rec2) in enumerate(
                itertools.izip(SeqIO.parse(misc.gzip_friendly_open(fpath1), 'fastq'),
                               SeqIO.parse(misc.gzip_friendly_open(fpath2), 'fastq'))
        ):
            total += 1
            seq = classify_seq(rec1, rec2)
            if seq:
                read_names_given_seq[seq].append(str(rec1.id))
            else:
                discarded += 1
        found = total - discarded
        log.info('    Found {} of {} ({:.1f}%)'.format(found, total, 100 * found / float(total)))

    #--------------------------------------------------------------------------------
    # Output results
    #--------------------------------------------------------------------------------
    out_fpath = '{}_read_names_by_seq.txt'.format(arguments.out_prefix)
    with open(out_fpath, 'w') as out:
        for seq, read_names in sorted(read_names_given_seq.items()):
            out.write('{}\t{}\n'.format(seq, '\t'.join(read_names)))


def make_read_names_by_sample(arguments):
    out_fpath = '{}_read_names_by_sample.txt'.format(arguments.out_prefix)
    log.info('Making {}'.format(out_fpath))
    with open(out_fpath, 'w') as out:
        for sample_dir in arguments.sample_dirs:
            log.info('    {}'.format(sample_dir))
            out.write('>{}\n'.format(sample_dir))
            fastq_fpaths = [fpath for fpath in glob.glob(sample_dir, '*')]
            pe_fpaths, se_fpaths = misc.find_paired_and_unpaired_files_from_fpaths(fastq_fpaths)
            assert not se_fpaths, 'All fastq files should be paired:\n{}'.format('\n'.join(se_fpaths))
            for fpath1, fpath2 in pe_fpaths:
                for rec in SeqIO.parse(misc.gzip_friendly_open(fpath1), 'fastq'):
                    out.write('{}\n'.format(str(rec.id)))
    log.info('Finished making {}'.format(out_fpath))


def setup_run(arguments):
    stitch_reads(arguments)
    make_read_names_by_sample(arguments)
