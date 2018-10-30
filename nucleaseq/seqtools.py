import os
import string
import logging
from collections import Counter, defaultdict

log = logging.getLogger(__name__)

bases = 'ACGT'

dna_complements = string.maketrans('acgtnACGTN', 'tgcanTGCAN')
def dna_rev_comp(dna_string):
    return dna_string.translate(dna_complements)[::-1]


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
            read_names_given_sample[sample] = line.strip()
    return read_names_given_sample
