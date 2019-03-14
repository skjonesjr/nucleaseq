import numpy as np
from freebarcodes import seqtools


golden_ratio = (np.float128(1) + np.sqrt(np.float128(5)))/np.float128(2.0)
golden_remainder = golden_ratio % np.float128(1)

def float2int(f, k):
    return int(f * 2**(2*k))

def float2dna(f, k):
    return seqtools.num2dna(float2int(f, k), k)


def passes_filters(seq, GC_max, bad_substrs=['CC', 'GG', 'AAA', 'TTT']):
    GC_count = seq.count('C') + seq.count('G')
    if GC_count > GC_max or len(seq) - GC_count > GC_max:
        return False
    for substr in bad_substrs:
        if substr in seq:
            return False
    for i in range(len(seq) - 6):
        if seqtools.dna_rev_comp(seq[i:i+3]) in seq[i+3:]:
            return False
    return True

def golden_iterator(bc_len):
    GC_max = min(range(bc_len), key=lambda x: abs(float(x)/bc_len-0.6))
    
    val = 0
    for _ in xrange(10000000):
        val += golden_remainder
        val %= 1
        val_int = float2int(val, bc_len)
        val_seq = seqtools.num2dna(val_int, bc_len)
        
        if passes_filters(val_seq, GC_max):
            yield val_seq
