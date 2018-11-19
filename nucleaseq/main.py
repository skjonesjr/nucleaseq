"""
NucleaSeq

Usage:
  nucleaseq stitch <min_len> <max_len> <out_prefix> <fastq_files>...
  nucleaseq preprocess <uncut_or_cut> <targets_file> <target_name> <pamtarg_pos> <exploded_oligos_file> <max_primer_err> <max_bc_err> <read_names_by_seq_file> <read_names_by_sample_file> <out_prefix> <start> <inc> <large_inc> [--nprocs=<nprocs>] [-v | -vv | -vvv]

Options:
  -h --help     Show this screen.
  --version     Show version.

Commands:
  stitch            Find ML stitched sequence from Illumina PE reads
  preprocess        Preprocess reads to determine identity and structure.

"""
import logging
import os
from nucleaseq.constants import VERSION
from nucleaseq.config import CommandLineArguments
from nucleaseq.ml_stitch_pe_reads import stitch_reads
from nucleaseq.preprocess import preprocess
from docopt import docopt


def main(**kwargs):
    docopt_args = docopt(__doc__, version=VERSION)
    arguments = CommandLineArguments(docopt_args, os.getcwd())

    log = logging.getLogger()
    handler = logging.StreamHandler()
    formatter = logging.Formatter("%(asctime)s   %(message)s", "%Y-%m-%d %H:%M:%S")
    handler.setFormatter(formatter)
    log.addHandler(handler)
    log.setLevel(arguments.log_level)
    log.debug(docopt_args)

    commands = {
        'stitch': stitch_reads,
        'preprocess': preprocess,
    }

    commands[arguments.command](arguments)


if __name__ == '__main__':
    main()
