import random
import numpy as np
from seqtools import bases
from freebarcodes.seqtools import *
from collections import Counter
from copy import deepcopy


def generate_clean_random_eqmarg_seqs(nseq, seqlen):
    """
    Generates a set of random sequences where all row and column marginals
    are equal between the four bases (or within one if total not mult of 4).
    
        nseq :int:    Number of desired sequences
        seqlen :int:  Desired length of sequences
    """
    while True:
        seqs, nsteps = generate_random_eqmarg_seqs_with_no_triples(nseq, seqlen)
        if is_good_seqset(seqs):
            return seqs

def equalish_base_counts(n):
    base_counts = [int(n)/4]*4
    for i in range(n - sum(base_counts)):
        base_counts[i] += 1
    random.shuffle(base_counts)
    return base_counts

def block_matrix_with_equalish_base_counts(n, m):
    M = np.zeros((n, m), dtype=np.uint8)
    seq_base_counts = equalish_base_counts(m)
    col_base_counts = equalish_base_counts(n)
    for i in range(4):
        row_start = sum(col_base_counts[:i])
        row_end = row_start + col_base_counts[i]
        for j in range(4):
            col_start = sum(seq_base_counts[:j])
            col_end = col_start + seq_base_counts[j]
            M[row_start:row_end, col_start:col_end] = (i + j) % 4
    return M

def find_triplet_i_j(M):
    for i in range(M.shape[0]):
        for j in range(M.shape[1] - 2):
            if M[i, j] == M[i, j + 1] == M[i, j + 2]:
                return i, j + 1
    return False

def swap_creates_triple(tmp, mid_idx, b):
    tmp[mid_idx] = b
    for i in range(len(tmp) - 3):
        if np.array_equal(tmp[i:i+3], b * np.ones((3,))):
            return True
    return False

def constant_marginal_swap(M, i, j):
    tmp_j_idx = min(2, j)
    tmp_j_start = max(0, j-2)
    i2_start = random.randint(0, M.shape[0]-1)
    j2_start = random.randint(0, M.shape[1]-1)
    b_tl = M[i, j]
    for i2 in range(i2_start, i2_start + M.shape[0]):
        i2 %= M.shape[0]
        b_bl = M[i2, j]
        if b_tl == b_bl:
            continue
        tmp = deepcopy(M[i2, tmp_j_start:j+3])
        if swap_creates_triple(tmp, tmp_j_idx, b_tl):
            continue            
        for j2 in range(j2_start, j2_start + M.shape[1]):
            j2 %= M.shape[1]
            tmp_j2_idx = min(2, j)
            tmp_j2_start = max(0, j2-2)
            b_tr = M[i, j2]
            b_br = M[i2, j2]
            tmp_tr = deepcopy(M[i, tmp_j2_start:j2+3])
            tmp_br = deepcopy(M[i2, tmp_j2_start:j2+3])
            if (b_tl == b_br and b_tr == b_bl 
                    and not swap_creates_triple(tmp_tr, tmp_j2_idx, b_br)
                    and not swap_creates_triple(tmp_br, tmp_j2_idx, b_tr)):
                M[i, j] = b_tr
                M[i2, j] = b_tl
                M[i, j2] = b_br
                M[i2, j2] = b_bl
                return

def random_constant_marginal_swap(M):
    i = random.randint(0, M.shape[0]-1)
    j = random.randint(0, M.shape[1]-1)
    return constant_marginal_swap(M, i, j)

def fix_triple_by_constant_marginal_swap(M):
    try:
        i, j = find_triplet_i_j(M)
    except:
        return False
    constant_marginal_swap(M, i, j)
    return True

def fix_all_triples_by_constant_marginal_swap(M):
    nsteps = 0
    while fix_triple_by_constant_marginal_swap(M):
        nsteps += 1
    return nsteps
        
def is_good_seqset(seqs):
    nseq = len(seqs)
    seqlen = len(list(seqs)[0])
    if len(set(seqs)) != len(seqs):
        return False  # repeats
    for seq in seqs:
        if any(3*base in seq for base in bases):
            return False  # has triples
        cntr = Counter(seq)
        if min(cntr.values()) < int(seqlen)/4 or max(cntr.values()) > int(seqlen)/4 + 1:
            return False  # Non-uniform base distribution per seq
    for i in range(seqlen):
        cntr = Counter(seq[i] for seq in seqs)
        if min(cntr.values()) < int(nseq)/4 or max(cntr.values()) > int(nseq)/4 + 1:
            return False  # Non-uniform base distribution per position
    return True
    
def seqs_given_index_matrix(M):
    return [''.join([bases[M[i, j]] for j in range(M.shape[1])]) for i in range(M.shape[0])]

def generate_random_eqmarg_seqs_with_no_triples(nseq, seqlen):
    M = block_matrix_with_equalish_base_counts(nseq, seqlen)
    for _ in xrange(nseq * seqlen):
        random_constant_marginal_swap(M)
    nsteps = fix_all_triples_by_constant_marginal_swap(M)
    return seqs_given_index_matrix(M), nsteps
