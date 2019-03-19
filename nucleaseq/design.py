import random
import seqtools


def shuffled_equalish_bases(slen, bad_substrs):
    n_copies = int(slen + 3)/4
    l = list(bases * n_copies)[:slen]
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
                   target,
                   bad_substrs,
                   fudge_factor,
                   abs_cannonical_cut_sites):

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
