import logging
import editdistance
from freebarcodes import editmeasures
from collections import Counter
from general_sequence_tools import dna_rev_comp
from NucleaSeqOligo import NucleaSeqOligo


log = logging.getLogger(__name__)


def assert_len_1_and_return_element(my_set):
    assert len(my_set) == 1, len(my_set)
    return list(my_set)[0]


class OligosContainer(object):
    """
    A class for the convenient storage and processing of a library of NucleaSeq oligo sequences.
    """
    def __init__(self,
                 exploded_fpath,
                 perfect_target,
                 perfect_pamtarg_one_target_pos):

        # Load oligos
        self.oligos = [NucleaSeqOligo(*line.strip().split('\t')) for line in open(exploded_fpath)]
        self.oligos_set = set(self.oligos)
        self.perfect_target = perfect_target
        log.info('Loaded {:,d} oligos ({:,d} unique)'.format(len(self.oligos),
                                                            len(self.oligos_set)))
        log.info('Oligo lengths: {}'.format(Counter(map(len, self.oligos))))

        # Find useful subsets
        self.perfect_target_oligos = [oligo for oligo in self.oligos 
                                      if oligo._target == self.perfect_target]
        self.target_oligos = [oligo for oligo in self.oligos if oligo._buffer_left]
        self.non_target_oligos = [oligo for oligo in self.oligos if not oligo._buffer_left]
        log.info('{:,d} Perfect target oligos'.format(len(self.perfect_target_oligos)))
        log.info('{:,d} Target oligos'.format(len(self.target_oligos)))
        log.info('{:,d} Non-target oligos'.format(len(self.non_target_oligos)))

        # Find primers
        self.cr_left = assert_len_1_and_return_element(set(oligo._cr_left for oligo in self.oligos))
        self.cr_right_rc = assert_len_1_and_return_element(set(oligo._cr_right for oligo in self.oligos))
        self.cr_right = dna_rev_comp(self.cr_right_rc)
        log.info('Left/Right Primer Seqs: {} / {}'.format(self.cr_left, self.cr_right))
        log.info('Left/Right Primer Edit Distance: {}'.format(
            editdistance.eval(self.cr_left, self.cr_right))
        )

        # Check buffers seqs
        self.buffer_lefts = Counter(oligo._buffer_left for oligo in self.target_oligos)
        self.buffer_rights = Counter(oligo._buffer_right for oligo in self.target_oligos)
        log.info('Target oligo left buffers: {}'.format(Counter(self.buffer_lefts)))
        log.info('Target oligo right buffers: {}'.format(Counter(self.buffer_rights)))

        # Find seq subsets by buffer
        most_common_buffer_left = sorted(self.buffer_lefts,
                                         key=self.buffer_lefts.get,
                                         reverse=True)[0]
        most_common_buffer_right = sorted(self.buffer_rights,
                                          key=self.buffer_rights.get,
                                          reverse=True)[0]
        self.perfect_target_and_buffer_oligos = [
            oligo for oligo in self.perfect_target_oligos
            if oligo._buffer_left == most_common_buffer_left
            and oligo._buffer_right == most_common_buffer_right
        ]
        self.alt_buffer_oligos = [
            oligo for oligo in self.target_oligos
            if oligo._buffer_left != most_common_buffer_left
            or oligo._buffer_right != most_common_buffer_right
        ]
        log.info('Most common left/right buffers: {} / {}'.format(
            most_common_buffer_left, most_common_buffer_right
        ))
        log.info('{:,d} Perfect target and most common buffer oligos'.format(
            len(self.perfect_target_and_buffer_oligos)
        ))
        log.info('{:,d} Alternate buffer oligos'.format(len(self.alt_buffer_oligos)))

        # Check right_buffer_buffers 
        log.info('Target right buffer buffers: {}'.format(
            Counter(oligo._right_buffer_buffer for oligo in self.target_oligos)
        ))

        # Process barcodes
        self.left_barcodes = set(oligo._barcode_left for oligo in self.oligos)
        self.right_barcodes = set(dna_rev_comp(oligo._barcode_right) for oligo in self.oligos)
        self.barcodes = self.left_barcodes | self.right_barcodes
        report_str = 'Barcodes not present in other side\'s barcode list: Left {} / Right {}'.format(
                len(self.left_barcodes - self.right_barcodes),
                len(self.right_barcodes - self.left_barcodes),
            )
        if len(self.left_barcodes ^ self.right_barcodes) > 2:
            log.warn(report_str)
        else:
            log.info(report_str)

        # Get bc_len
        self.bc_len = assert_len_1_and_return_element(set(len(bc) for bc in self.barcodes))
        log.info('Barcode length: {}'.format(self.bc_len))

        # Make relevant dicts
        self.oligo_given_left_barcode = {oligo._barcode_left: oligo for oligo in self.oligos}
        self.oligo_given_right_barcode = {dna_rev_comp(oligo._barcode_right): oligo 
                                          for oligo in self.oligos}
        self.oligo_given_barcode_given_side = {'left': self.oligo_given_left_barcode, 
                                               'right': self.oligo_given_right_barcode}
        self.oligo_given_seq = {oligo.sequence: oligo for oligo in self.oligos}
        self.target_given_seq = {oligo.sequence: oligo._target for oligo in self.oligos}
        log.info('Made barcode and seq dicts.')

        # Update pamtarg coords 
        for oligo in self.target_oligos:
            oligo.add_pamtarg_coord_one_pos(perfect_target, perfect_pamtarg_one_target_pos)

        log.info('Pam-target-one target positions: {}'.format(
            Counter([oligo.pamtarg_one_target_pos for oligo in self.target_oligos]))
        )
