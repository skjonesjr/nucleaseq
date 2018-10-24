import logging
import os


class CommandLineArguments(object):
    """
    Wraps the raw arguments provided by docopt.
    """
    def __init__(self, arguments, current_directory):
        self._arguments = arguments
        self._current_directory = current_directory

    def _comma_delimited_arg(self, key):
        if self._arguments[key]:
            return self._arguments[key].split(',')
        return []

    @property
    def command(self):
        if self._arguments.get('preprocess'):
            return 'preprocess ' + self._arguments.get('<uncut_or_cut>')
        ## We have to do this weird loop to deal with the way docopt stores the command name
        #for possible_command in ('decode',
        #                         'generate',
        #                         'prune',
        #                         'concatenate'):
        #    if self._arguments.get(possible_command):
        #        return possible_command

    @property
    def targets_file(self):
        return os.path.expanduser(self._arguments['<targets_file>']) or None

    @property
    def exploded_oligos_file(self):
        return os.path.expanduser(self._arguments['<exploded_oligos_file>']) or None

    @property
    def read_names_by_seq_file(self):
        return os.path.expanduser(self._arguments['<read_names_by_seq_file>']) or None

    @property
    def read_names_by_sample_file(self):
        return os.path.expanduser(self._arguments['<read_names_by_sample_file>']) or None

    @property
    def target_name(self):
        return self._arguments['<target_name>']

    @property
    def pamtarg_pos(self):
        return int(self._arguments['<pamtarg_pos>'])

    @property
    def max_primer_err(self):
        return int(self._arguments['<max_primer_err>'])

    @property
    def max_bc_err(self):
        return int(self._arguments['<max_bc_err>'])

    @property
    def nprocs(self):
        if self._arguments['--nprocs']:
            return int(self._arguments['--nprocs'])
        return 1

    @property
    def start(self):
        return int(self._arguments['<start>'] or 0)

    @property
    def inc(self):
        return int(self._arguments['<inc>'] or 100000)

    @property
    def large_inc(self):
        return int(self._arguments['<large_inc>'] or 1000000)

    @property
    def log_level(self):
        log_level = {0: logging.ERROR,
                     1: logging.WARN,
                     2: logging.INFO,
                     3: logging.DEBUG}
        # default to silent if the user supplies no verbosity setting
        return log_level.get(self._arguments['-v'], logging.ERROR)
