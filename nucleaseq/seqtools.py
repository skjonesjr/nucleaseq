import os
import string
import logging
from collections import defaultdict

log = logging.getLogger(__name__)

bases = 'ACGT'

dna_complements = string.maketrans('acgtnACGTN', 'tgcanTGCAN')
def dna_rev_comp(dna_string):
    return dna_string.translate(dna_complements)[::-1]
def dna_fwd_comp(dna_string):
    return dna_string.translate(dna_complements)


def load_deduped_read_names_given_seq(fpath):
    log.info('Loading read_names_given_seq file {}'.format(fpath))
    read_names_given_seq = defaultdict(set)
    for line in open(fpath):
        words = line.strip().split()
        seq = words[0]
        read_names = words[1:]
        read_names_given_seq[seq].update(read_names)
    read_names_given_seq = dict(read_names_given_seq)
    log.info('Raw read_names_given_seq: {:,d} seqs, {:,d} reads'.format(
        len(read_names_given_seq),
        sum(len(v) for v in read_names_given_seq.values())
    ))

    # De-duplicating reverse-comp seqs
    for seq, rns in read_names_given_seq.items():
        if seq not in read_names_given_seq:
            continue
        seq_rc = dna_rev_comp(seq)
        if seq_rc in read_names_given_seq:
            read_names_given_seq[seq].update(read_names_given_seq[seq_rc])
            del read_names_given_seq[seq_rc]
    log.info('Deduped read_names_given_seq: {:,d} seqs, {:,d} reads'.format(
        len(read_names_given_seq),
        sum(len(v) for v in read_names_given_seq.values())
    ))
    return read_names_given_seq


def load_read_name_seq_items(fpath):
    suffix = '.dedup_sort.txt'
    if fpath.endswith(suffix):
        out_fpath = fpath
    else:
        out_fpath = fpath + suffix

    if os.path.exists(out_fpath):
        log.info('Loading deduped and sorted read_names_given_seq: {}'.format(out_fpath))
        read_name_seq_items = []
        for line in open(out_fpath):
            words = line.strip().split()
            read_name_seq_items.append((words[0], words[1:]))
    else:
        read_names_given_seq = load_deduped_read_names_given_seq(fpath)
        # Saving as items and sorting
        log.info('Sorting read_name_seq_items')
        read_name_seq_items = read_names_given_seq.items()
        read_name_seq_items.sort()
        with open(out_fpath, 'w') as out:
            for seq, rns in read_name_seq_items:
                out.write('{}\t{}\n'.format(seq, '\t'.join(rns)))

    log.info('Deduped sorted read_name_seq_items: {:,d} seqs, {:,d} reads'.format(
        len(read_name_seq_items),
        sum(len(tup[1]) for tup in read_name_seq_items)
    ))
    return read_name_seq_items


def load_sample_given_read_name(fpath):
    log.info('Loading read_names_by_sample file {}'.format(fpath))
    sample_given_read_name = {}
    for line in open(fpath):
        if line.startswith('>'):
            sample = line.strip()[1:]
        else:
            sample_given_read_name[line.strip()] = sample
    return sample_given_read_name


def load_read_names_given_sample(fpath):
    log.info('Loading read_names_given_sample from {}'.format(fpath))
    read_names_given_sample = defaultdict(list)
    for line in open(fpath):
        if line.startswith('>'):
            sample = line.strip()[1:]
        else:
            read_names_given_sample[sample].append(line.strip())
    return read_names_given_sample


def iterate_single_mismatches(perfect_target):
    for si in range(len(perfect_target)):
        for bi, base_i in enumerate(bases.replace(perfect_target[si], '')):
            seq = perfect_target[:si] + base_i + perfect_target[si+1:]
            yield si, bi, base_i, seq


def iterate_double_mismatches(perfect_target):
    for sj in range(len(perfect_target)):
        for si in range(sj):
            for bi, base_i in enumerate(bases.replace(perfect_target[si], '')):
                for bj, base_j in enumerate(bases.replace(perfect_target[sj], '')):
                    seq = perfect_target[:si] + base_i + perfect_target[si+1:sj] + base_j + perfect_target[sj+1:]
                    yield si, sj, bi, bj, base_i, base_j, seq


def iterate_single_deletions(perfect_target):
    for si in range(len(perfect_target)):
        seq = perfect_target[:si] + perfect_target[si+1:]
        yield si, seq


def iterate_double_deletions(perfect_target):
    for sj in range(len(perfect_target)):
        for si in range(sj):
            seq = perfect_target[:si] + perfect_target[si+1:sj] + perfect_target[sj+1:]
            yield si, sj, seq


def iterate_single_insertions(perfect_target):
    for si in range(len(perfect_target) + 1):
        for bi, base_i in enumerate(bases):
            seq = perfect_target[:si] + base_i + perfect_target[si:]
            yield si, bi, seq


def iterate_double_insertions(perfect_target):
    for sj in range(len(perfect_target) + 1):
        for si in range(sj):
            for bi, base_i in enumerate(bases):
                for bj, base_j in enumerate(bases):
                    seq = perfect_target[:si] + base_i + perfect_target[si:sj] + base_j + perfect_target[sj:]
                    yield si, sj, bi, bj, seq


def iterate_complement_stretches(perfect_target):
    for end in range(len(perfect_target)+1):
        for start in range(end):
            seq = perfect_target[:start] + dna_fwd_comp(perfect_target[start:end]) + perfect_target[end:]
            yield start, end, seq


def double_mismatch_matrix(val_given_seq, perfect_target):
    s = 3 * len(perfect_target)
    M = np.zeros((s, s))
    M[:, :] = None
    for si, sj, bi, bj, base_i, base_j, seq in iterate_double_mismatches(perfect_target):
        M[3 * sj + bj, 3 * si + bi] = val_given_seq.get(seq, None)
    for si, bi, base_i, seq in iterate_single_mismatches(perfect_target):
        idx = 3*si + bi
        M[idx, idx] = val_given_seq.get(seq, None)
    return M


def double_deletion_matrix(val_given_seq, perfect_target):
    s = len(perfect_target)
    M = np.zeros((s, s))
    M[:, :] = None
    for si, sj, seq in iterate_double_deletions(perfect_target):
        M[sj, si] = val_given_seq.get(seq, None)
    for si, seq in iterate_single_deletions(perfect_target):
        M[si, si] = val_given_seq.get(seq, None)
    return M


def double_insertion_matrix(val_given_seq, perfect_target):
    s = 4 * (len(perfect_target) + 1)
    M = np.zeros((s, s))
    M[:, :] = None
    for si, sj, bi, bj, seq in iterate_double_insertions(perfect_target):
        M[4 * sj + bj, 4 * si + bi] = val_given_seq.get(seq, None)
    for si, bi, seq in iterate_single_insertions(perfect_target):
        idx = 4*si + bi
        M[idx, idx] = val_given_seq.get(seq, None)
    return M


def complement_stretches_matrix(val_given_seq, perfect_target):
    s = len(perfect_target)
    M = np.zeros((s, s))
    M[:, :] = None
    for start, end, seq in iterate_complement_stretches(perfect_target):
        M[end - 1, start] = val_given_seq.get(seq, None)
    return M
