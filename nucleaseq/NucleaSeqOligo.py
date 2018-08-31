from Bio import pairwise2
from align import aligner
from align.matrix import DNAFULL

class NucleaSeqOligo(object):
    """
    A class for the convenient storage and processing of one NucleaSeq oligo sequence.
    """
    self.pieces_names = [
        'primer_left',
        'barcode_left',
        'buffer_left',
        'target',
        'buffer_right',
        'right_buffer_buffer',
        'barcode_right',
        'primer_right',
    ]

    def __init__(self,
                 cr_left,
                 barcode_left,
                 buffer_left,
                 target,
                 buffer_right,
                 right_buffer_buffer,
                 barcode_right,
                 cr_right):
        self._cr_left = cr_left
        self._barcode_left = barcode_left
        self._buffer_left = buffer_left
        self._target = target
        self._buffer_right = buffer_right
        self._right_buffer_buffer = right_buffer_buffer
        self._barcode_right = barcode_right
        self._cr_right = cr_right

    @property
    def pieces(self):
        return [self._cr_left,
                self._barcode_left,
                self._buffer_left,
                self._target,
                self._buffer_right,
                self._right_buffer_buffer,
                self._barcode_right,
                self._cr_right]

    @property
    def pieces_and_names(self):
        return zip(self.pieces, self.pieces_names)

    @property
    def sequence(self):
        return ''.join(self.pieces)
    @property
    def gc_content(self):
        sequence = self.sequence
        g = sequence.count('G')
        c = sequence.count('C')
        return float(g + c) / float(len(sequence))
    
    def __hash__(self):
        return hash(self.sequence)
    
    def __eq__(self, other):
        return self.sequence == other.sequence
    
    def __len__(self):
        return len(self.sequence)
    
    def add_right_buffer_buffer(self, right_buffer_buffer, desired_oligo_len):
        rbb_len = desired_oligo_len - (len(self) - len(self._right_buffer_buffer))
        self._right_buffer_buffer = right_buffer_buffer[:rbb_len]
        
    def tab_delimited_str(self):
        return '\t'.join(self.pieces)
    
    def add_pamtarg_coord_one_pos(self, perfect_target, perfect_pamtarg_one_target_pos):
        # Find position in the current oligo corresponding to "1" in pam-target coordinates. That
        # is, where the first base after the pam-target boundary is "1", the first before the
        # boundary is "-1", etc.
        alignments = pairwise2.align.globalms(perfect_target, self._target, 2, -1, -1, -0.9)
        ref_align, observed_align = alignments[0][:2]

        ref_idx = 0
        obs_idx = 0
        for rc, oc in zip(ref_align, observed_align):
            if ref_idx == perfect_pamtarg_one_target_pos:
                self.pamtarg_one_target_pos = obs_idx
                self.pamtarg_one_oligo_pos = obs_idx + sum(len(piece) for piece in self.pieces[:3])
                return
            if rc != '-':
                ref_idx += 1
            if oc != '-':
                obs_idx += 1

        raise ValueError('Could not find one in pamtarg coords:{}\n{}\n{}\n'.format(
            perfect_pamtarg_one_target_pos,
            ref_align,
            observed_align,
        ))

    def convert_target_pos_to_pamtarg_coord(self, pos):
        pamtarg_coord = pos - self.pamtarg_one_target_pos
        if pamtarg_coord >= 0:
            pamtarg_coord += 1
        return pamtarg_coord

    def convert_oligo_pos_to_pamtarg_coord(self, pos):
        pamtarg_coord = pos - self.pamtarg_one_oligo_pos
        if pamtarg_coord >= 0:
            pamtarg_coord += 1
        return pamtarg_coord
    
    def convert_pamtarg_coord_to_target_pos(self, pamtarg_coord):
        pos = pamtarg_coord + self.pamtarg_one_target_pos
        if pamtarg_coord >= 0:
            pos -= 1
        return pos

    def convert_pamtarg_coord_to_oligo_pos(self, pamtarg_coord):
        pos = pamtarg_coord + self.pamtarg_one_oligo_pos
        if pamtarg_coord >= 0:
            pos -= 1
        return pos

    def prefix_to_pamtarg_coord(self, pamtarg_coord):
        pos = self.convert_pamtarg_coord_to_oligo_pos(pamtarg_coord)
        return self.sequence[:pos]

    def suffix_to_pamtarg_coord(self, pamtarg_coord):
        pos = self.convert_pamtarg_coord_to_oligo_pos(pamtarg_coord)
        return self.sequence[pos:]
    
    def align_and_return_as_pieces(self, observed_seq, align_method='global'):
        # For cut products, best method is global_cfe (Cost Free Ends)
        if align_method == 'global':
            alignments = pairwise2.align.globalms(self.sequence, observed_seq, 2, -1, -1, -0.9)
            ref_align, observed_align = alignments[0][:2]  
        elif align_method == 'global_cfe':
            ref_align, observed_align = aligner(self.sequence, observed_seq, method='global_cfe',
                                                matrix=DNAFULL, gap_open=-10, gap_extend=-1)[0][:2]
        else:
            raise ValueError(
                'Only supported alignment methods are global / global_cfe: {}'.format(align_method)
            )

        ref_idx = 0
        piece_idx = 0
        starts = [0]
        for i, (rc, oc) in enumerate(zip(ref_align, observed_align)):
            if rc != '-':
                while ref_idx == sum(len(p) for p in self.pieces[:piece_idx + 1]): # while takes care of length-zero pieces
                    starts.append(i)
                    piece_idx += 1
                ref_idx += 1
        assert len(starts) == len(self.pieces), (ref_align, observed_align,
                                                       self.pieces, starts)

        ends = starts[1:] + [len(ref_align)]
        return [(ref_align[start:end], observed_align[start:end])
                for start, end in zip(starts, ends)]
