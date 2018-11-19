import os
import gzip
import re


def gzip_friendly_open(fname, mode='r'):
    """
    gzip_friendly_open returns a file handle appropriate for the file, whether gzipped or not.
    """
    if os.path.splitext(fname)[-1] == '.gz':
        return gzip.open(fname, mode)
    else:
        return open(fname, mode)


def find_paired_and_unpaired_files_from_fpaths(fq_fpaths, skip_I_files=True):
    se_fpaths = []
    pe_fpaths = []
    fq_fpaths.sort()
    i = 0
    while i < len(fq_fpaths):
        if skip_I_files and '_I1_' in fq_fpaths[i] or '_I2_' in fq_fpaths[i]:
            i += 1
            continue
        if i + 1 == len(fq_fpaths):
            se_fpaths.append(fq_fpaths[i])
            break
        m1 = re.match('^(.*)_R1(.*)$', fq_fpaths[i])
        m2 = re.match('^(.*)_R2(.*)$', fq_fpaths[i+1])
        if m1 and m2 and m1.group(1) == m2.group(1) and m1.group(2) == m2.group(2):
            pe_fpaths.append((fq_fpaths[i], fq_fpaths[i+1]))
            i += 2
        else:
            se_fpaths.append(fq_fpaths[i])
            i += 1
    return pe_fpaths, se_fpaths
