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
