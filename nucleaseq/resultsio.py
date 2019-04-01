import os
from collections import namedtuple


def write_cleavage_rates_file(run_label,
                              oligo_container,
                              vk_y0_given_oligo_seq,
                              bootstrap_vk_y0s,
                              dpath='.'):
    cut_fpath = os.path.join('{}_cleavage_rate_and_y0.txt'.format(run_label))
    with open(cut_fpath, 'w') as out:
        oligo = oligo_container.oligo_given_seq[vk_y0_given_oligo_seq.keys()[0]]
        out.write('\t'.join(['\t'.join(oligo.pieces_names),
                             'cleavage_rate',
                             'log10_cleavage_rate_err',
                             'cleavage_rate_5th_pctl',
                             'cleavage_rate_95th_pctl',
                             'y0',
                             'y0_err']) + '\n')
        for seq, (vk, y0) in sorted(vk_y0_given_oligo_seq.items()):
            other_vks = [d_vk_y0[seq][0] for d_vk_y0 in bootstrap_vk_y0s]
            other_y0s = [d_vk_y0[seq][1] for d_vk_y0 in bootstrap_vk_y0s]
            log_vk_err = np.std([np.log10(vk) for vk in other_vks])
            pctl5 = np.percentile(other_vks, 5)
            pctl_95 = np.percentile(other_vks, 95)
            y0_err = np.std(other_y0s)
            oligo = oligo_container.oligo_given_seq[seq]
            out.write('\t'.join(map(str, ['\t'.join(oligo.pieces), vk, log_vk_err, pctl5, pctl_95, y0, y0_err])) + '\n')


def read_cleavage_rates_file(fpath):
    records = []
    with open(fpath) as f:
        headers = next(f).strip().split()
        cr_idx = headers.index('cleavage_rate')
        OligoRecord = namedtuple('OligoRecord', 'oligo_seq ' + ' '.join(headers))
        for line in f:
            words = line.strip().split('\t')
            vals = words[:cr_idx] + [float(w) for w in words[cr_idx:]]
            oligo_seq = ''.join([words[i] for i in range(cr_idx)])
            records.append(OligoRecord(oligo_seq, *vals))
    return records
