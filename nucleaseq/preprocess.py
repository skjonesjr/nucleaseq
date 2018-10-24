import os
from subprocess import check_call
from seqtools import load_deduped_read_names_given_seq
import logging

log = logging.getLogger(__name__)

THIS_DIR = os.path.dirname(os.path.abspath(__file__))

def preprocess(arguments):
    read_name_items = load_deduped_read_names_given_seq(read_names_by_seq_file)
    last_end = len(read_name_items)

    if arguments.cut_or_uncut == 'cut':
        preprocess_script = os.path.join(THIS_DIR, 'preprocess_cut.py')
    if arguments.cut_or_uncut == 'uncut':
        preprocess_script = os.path.join(THIS_DIR, 'preprocess_uncut.py')

    for start in range(0, last_end, arguments.large_inc):
        cmd = [
            'python',
            preprocess_script,
            arguments.targets_file,
            arguments.target_name,
            arguments.pamtarg_pos,
            arguments.exploded_oligos_file,
            arguments.max_primer_err,
            arguments.max_bc_err,
            arguments.read_names_by_seq_file,
            arguments.read_names_by_sample_file,
            arguments.start,
            arguments.inc,
            arguments.start + arguments.large_inc,
            arguments.nprocs,
        ]
        cmd = map(str, cmd)
        check_call(cmd)
