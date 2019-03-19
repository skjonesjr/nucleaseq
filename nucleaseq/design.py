import sys
import random
import seqtools
import editdistance
from golden_iterator import parallel_golden_iterator
from multiprocessing import Pool


def shuffled_equalish_bases(slen, bad_substrs):
    n_copies = int(slen + 3)/4
    l = list(seqtools.bases * n_copies)[:slen]
    while True:
        random.shuffle(l)
        s = ''.join(l)[:slen]
        if not any(bad_sub in s for bad_sub in bad_substrs):
            return s


def no_runs(s):
    for c1, c2 in zip(s, s[1:]):
        if c1 == c2:
            return False
    return True


def get_buffer(slen, side, target, bad_substrs):
    assert side in ['left', 'right'], side
    while True:
        buf = shuffled_equalish_bases(slen, bad_substrs)
        if side == 'left' and buf[-1] != target[0] and buf[-1] not in 'CG' and no_runs(buf):
            return buf
        if side == 'right' and buf[0] != target[-1] and buf[0] not in 'CG' and no_runs(buf):
            return buf


def update_buffers(complete_sequences,
                   primer_len,
                   min_buffer_len,
                   target,
                   bad_substrs,
                   fudge_factor,
                   abs_cannonical_cut_sites,
                   pamtarg_coord_one_pos):

    # Left buffer
    left_len = max(min_buffer_len, max(primer_len - (ccs - fudge_factor) + 1
                                       for ccs in abs_cannonical_cut_sites))
    left_buffer = get_buffer(left_len, 'left', target, bad_substrs)
    for oligo in complete_sequences:
        if oligo._buffer_left:
            oligo._buffer_left = left_buffer

    # Right buffer
    right_len = max(min_buffer_len, max(primer_len - (len(target) - (ccs + fudge_factor)) + 1 
                                        for ccs in abs_cannonical_cut_sites))
    right_buffer = get_buffer(right_len, 'right', target, bad_substrs)
    for oligo in complete_sequences:
        if oligo._buffer_right:
            oligo._buffer_right = right_buffer
        
    # Right buffer buffer
    for oligo in complete_sequences:
        oligo._right_buffer_buffer = ''
    raw_oligo_lens = set(map(len, complete_sequences))
    if len(set(raw_oligo_lens)) == 1:
        return complete_sequences
    right_buffer_buffer = shuffled_equalish_bases(
        max(raw_oligo_lens) - min(raw_oligo_lens), 
        bad_substrs
    )
    desired_oligo_len = max(raw_oligo_lens)
    for oligo in complete_sequences:
        oligo.add_right_buffer_buffer(right_buffer_buffer, desired_oligo_len)
    oligo_lens = set(map(len, complete_sequences))
    assert len(oligo_lens) == 1, oligo_lens

    # Reset pamtarg one pos
    for oligo in complete_sequences:
        oligo.add_pamtarg_coord_one_pos(target, pamtarg_coord_one_pos)    
    
    return complete_sequences

def get_cut_prefixes(complete_sequences, cannonical_cut_sites, fudge_factor):
    cut_prefixes = set()
    for site in cannonical_cut_sites:
        for pamtarg_coord in range(site - fudge_factor, site + fudge_factor + 1 + 1):
            for oligo in complete_sequences:
                cut_prefixes.add(seqtools.dna_rev_comp(oligo.prefix_to_pamtarg_coord(pamtarg_coord)))
                cut_prefixes.add(oligo.suffix_to_pamtarg_coord(pamtarg_coord))
    return cut_prefixes


def is_desired_distance_from_all(seq, prefixes, dist):
    for prefix in prefixes:
        if editdistance.eval(seq, prefix) < dist:
            return False
    return True


def seq_if_potential_good_prefix(params):
    seq, primer_len_prefixes, cut_prefix_min_dist = params
    if is_desired_distance_from_all(seq, primer_len_prefixes, cut_prefix_min_dist):
        return seq


def find_good_prefixes(complete_sequences,
                       primer_len,
                       bad_substrs,
                       cannonical_cut_sites,
                       fudge_factor,
                       nprocs,
                       primer_max_err=2,
                       chunk_size=500000):
    cut_prefixes = get_cut_prefixes(complete_sequences, cannonical_cut_sites, fudge_factor)
    primer_len_prefixes = set([prefix[:primer_len] for prefix in cut_prefixes])
    primer_iter = parallel_golden_iterator(primer_len, nprocs, bad_substrs, max_tries=float('inf'))
    cut_prefix_min_dist = 4 * primer_max_err + 1
    good_prefix_min_dist = cut_prefix_min_dist + 2

    def add_to_good_prefixes(good_prefixes):
        next_seqs = [(next(primer_iter), primer_len_prefixes, cut_prefix_min_dist)
                     for _ in range(chunk_size)]
        sys.stdout.write('.')
        sys.stdout.flush()

        pl = Pool(nprocs)
        res = pl.map(seq_if_potential_good_prefix, next_seqs)
        pl.close()
        potential_good_prefixes = [seq for seq in res if seq]

        for seq in potential_good_prefixes:
            if is_desired_distance_from_all(seq, good_prefixes, good_prefix_min_dist):
                sys.stdout.write('*')
                sys.stdout.flush()
                good_prefixes.append(seq)
        return good_prefixes

    good_prefixes = []
    good_prefixes = add_to_good_prefixes(good_prefixes)
    i = 0
    while len(good_prefixes) == 1 and i < 10: # Accept 2+, give up if 0 or 10+ tries
        i += 1
        sys.stdout.write('.')
        sys.stdout.flush()
        good_prefixes = add_to_good_prefixes(good_prefixes)
    return good_prefixes
