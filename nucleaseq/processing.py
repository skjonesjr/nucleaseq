import os
import re
import numpy as np
import pandas as pd
import logging

log = logging.getLogger()
if not log.handlers:
    handler = logging.StreamHandler()
    formatter = logging.Formatter("%(asctime)s   %(message)s", "%Y-%m-%d %H:%M:%S")
    handler.setFormatter(formatter)
    log.addHandler(handler)
    log.setLevel(logging.INFO)


def load_and_concat_dfs(dpath, label):
    name_re = re.compile('{}.(\d+)-(\d+).pkl'.format(label))
    fnames = [fname for fname in os.listdir(dpath) if name_re.search(fname)]
    fnames.sort(key = lambda s: int(name_re.search(s).group(1)))
    
    m = name_re.search(fnames[0])
    start1, end1 = map(int, (m.group(1), m.group(2)))
    m = name_re.search(fnames[-1])
    start2, end2 = map(int, (m.group(1), m.group(2)))
    log.info('{}: {:,d}-{:,d}  ({:,d} files)'.format(label, start1, end2, len(fnames)))
    last_end = end2
    
    gap_total = 0
    for fname1, fname2 in zip(fnames, fnames[1:]):
        m = name_re.search(fname1)
        start1, end1 = map(int, (m.group(1), m.group(2)))
        m = name_re.search(fname2)
        start2, end2 = map(int, (m.group(1), m.group(2)))
        gap = start2 - end1 - 1
        if not gap == 0:
            log.info('Missing file between {} {} ({:,d} missing)'.format(fname1, fname2, gap))
            gap_total += gap
    log.info('Total missing: {:,d}'.format(gap_total))
            
    dfs = [pd.read_pickle(os.path.join(dpath, fname)) for fname in fnames]
    return pd.concat(dfs)


def calculate_y0s_and_vks(target_tab, label, bootstrap=False, verbose=False):
    # TODO clean this up
    assert all(target_tab.index == target_tab_multinom.index[:-1])
    cut_success = False
    while not cut_success:
        if bootstrap:
            if verbose: 
                print 'Calculating bootstrap read counts...'

            for sample in sample_names:
                target_tab[sample] = np.random.multinomial(target_tab_totals[sample], target_tab_multinom[sample])[:-1]

        if verbose: 
            print 'Normalizing...'

        norm_target_tab = target_tab.div(normalization_factors, axis='columns')
        norm_target_control_max = norm_target_tab.max(axis='index')[0]
        renorm_target_tab = norm_target_tab.div(norm_target_control_max)

        norm_perfect_target_tab = norm_target_tab.loc[[oligo.sequence for oligo in oc.perfect_target_oligos]]
        norm_perfect_target_control_max = norm_perfect_target_tab.max(axis='index')[0]
        renorm_perfect_target_tab = norm_perfect_target_tab.div(norm_perfect_target_control_max)

        nontarget_tab = pd.crosstab(non_target_data.oligo, non_target_data['sample'])
        norm_nontarget_tab = nontarget_tab.div(normalization_factors, axis='columns')
        norm_nontarget_control_max = norm_nontarget_tab.max(axis='index')[0]
        renorm_nontarget_tab = norm_nontarget_tab.div(norm_nontarget_control_max)

        non_target_medians = list(norm_nontarget_tab.median())
        Z = np.array(non_target_medians) / non_target_medians[0]
        med_norm_target_tab = norm_target_tab.div(Z, axis='columns')
        med_norm_perfect_target_tab = norm_perfect_target_tab.div(Z, axis='columns')
        med_norm_nontarget_tab = norm_nontarget_tab.div(Z, axis='columns')

        if verbose:
            print 'Reorganizing...'

        cut_target_tab = target_tab[cut_samples]

        med_norm_cut_target_tab = med_norm_target_tab[cut_samples]

        if verbose:
            print 'Calculating exp floors...'
        all_final_fracs = []
        exp_floors = []
        samples = cut_samples  # TODO 
        df = med_norm_perfect_target_tab.loc[:, samples]
        final_fracs = []
        for index, row in df.iterrows():
            y = np.array(row)
            y /= y[0]
            final_fracs.append(y[-1])
        plot_cdf(ax, final_fracs, label=label)
        all_final_fracs.extend(final_fracs)
        exp_floors.append(np.median(final_fracs))

        def exponential(x, y0, vk):
            return y0 * ((1 - exp_floors[0]) * np.exp(-vk * x) + exp_floors[0])

        def find_vk_0(seq_read_counts):
            halfway = 0.5 * seq_read_counts[0]
            for tA, tB, countA, countB in zip(times, times[1:], seq_read_counts, seq_read_counts[1:]):
                if countB <= halfway:
                    assert countA > halfway, counts
                    t_half = tA + (tB - tA) * (halfway - countA)/(countB - countA)
                    return np.log(2)/t_half
            return np.log(2)/times[-1]

        def curve_fit_cut_data(seq_read_counts):
            seq_read_counts = np.array(seq_read_counts) + 0.1
            y0_0 = seq_read_counts[0]
            vk_0 = find_vk_0(seq_read_counts)
        #    print y0_0, np.log(2)/vk_0
            return curve_fit(exponential, times, seq_read_counts, p0=[y0_0, vk_0], bounds=((0, 0), (1000, 1000)),
                             method='trf', max_nfev=10000)

        if verbose:
            print 'First fitting...'

        vks, y0s = [], []
        covs = []
        cut_success = True
        for i in range(len(target_tab)):
            try:
                test_counts = np.array(list(med_norm_cut_target_tab.iloc[i]))
                popt, pcov = curve_fit_cut_data(test_counts)
                y0s.append(popt[0])    
                vks.append(popt[1])
                covs.append(pcov)
            except:
                cut_success = False
                sys.stdout.write('^')
                break
        if not cut_success:
            continue

        vk_y0_given_oligo_seq = {}
        kkm_given_target = {}
        for oligo_seq, y0, vk in zip(med_norm_cut_target_tab.index, y0s, vks):
            vk_y0_given_oligo_seq[oligo_seq] = (vk, y0)

        if verbose:
            print 'Done'
        return vk_y0_given_oligo_seq, ampl_vk_y0_given_oligo_seq
