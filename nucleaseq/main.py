"""
NucleaSeq

Usage:
  nucleaseq setup_run <min_seq_len> <max_seq_len> <out_prefix> <sample_dirs>... [-v | -vv | -vvv]
  nucleaseq preprocess <uncut_or_cut> <targets_file> <target_name> <pamtarg_pos> <exploded_oligos_file> <max_primer_err> <max_bc_err> <read_names_by_seq_file> <read_names_by_sample_file> <out_prefix> <start> <inc> <large_inc> [--nprocs=<nprocs>] [-v | -vv | -vvv]

Commands:
  setup_run         Makes read_names_by_seq and read_names_by_sample for a given run.
  preprocess        Preprocess reads to determine identity and structure.
  
Parameters:
  min_seq_len                    Minimum sequence length to analyse. 
                                 Set to <= length of your shortest library member after cleavage but before adapter ligation.
  max_seq_len                    Maximum sequence length to analyse.
                                 Set to >= length of your largest library member before adapter ligation.
  out_prefix                     Prefix for output files for specific run
  sample_dirs                    Directories with sample fastq files, one subdirectory per time point. The directory path should end /* to encompass all sample subdirectories.
                                 Fastq files within the subdirectories may have .txt or .fastq extensions, and can be compressed (.gz).
                                 Paired reads in each subdirectory should have fastq filenames containing "_R1" and "_R2".
  uncut_or_cut                   "cut" or "uncut", selecting which reads to process
  targets_file                   Yaml file with targets and target names (see targets.yml)
  target_name                    Target name for current experiment
  pamtarg_pos                    Position in target defined as 1. Previous is -1.
  exploded_oligos_file           Exploded oligos file from library design
  max_primer_err                 Maximum allowed primer error
  max_bc_err                     Maximum allowed barcode error
  read_names_by_seq_file         Path to file generated during setup_run
  read_names_by_sample_file      Path to file generated during setup_run
  start                          Read index to start analysis (typically zero)
  inc                            Chunk size for data output (typically 100,000)
  large_inc                      Processing increment size (typically 2,000,000)
  nprocs                         Number of processors for parallelization

Options:
  -h --help     Show this screen.
  --version     Show version.

"""
import logging
import os
from nucleaseq.constants import VERSION
from nucleaseq.config import CommandLineArguments
from nucleaseq.setup_run import setup_run
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
        'setup_run': setup_run,
        'preprocess': preprocess,
    }

    commands[arguments.command](arguments)


if __name__ == '__main__':
    main()
