'''
#class _BedBuffer(DelimitedBuffer):
#    dataclass = Interval
#
#    def get_intervals(self):
#        self.validate_if_not()
#        chromosomes = self.get_text(0, fixed_length=False)
#        positions = self.get_integers(cols=[1, 2])
#        return Interval(chromosomes, positions[..., 0], positions[..., 1])
#
#    @classmethod
#    def _from_data(cls, data):
#        start_lens = np.log10(data.start).astype(int)+1
#        end_lens = np.log10(data.end).astype(int)+1
#        chromosome_lens = data.chromosome.shape.lengths
#        line_lengths = chromosome_lens + 1 + start_lens + 1 + end_lens + 1
#        line_ends = np.cumsum(line_lengths)
#        buf = np.empty(line_ends[-1], dtype=np.uint8)
#        lines = RaggedArray(buf, line_lengths)
#        obj = cls(buf, line_ends-1)
#        obj._move_2d_array_to_intervals(cls._move_ints_to_digit_array(data.end, np.max(end_lens)),
#                                        line_ends-1-end_lens, line_ends-1)
#
#        obj._move_2d_array_to_intervals(cls._move_ints_to_digit_array(data.start, np.max(start_lens)),
#                                        line_ends-2-end_lens-start_lens, line_ends-2-end_lens)
#
#        indices, _ = RaggedView(lines.shape.starts, chromosome_lens).get_flat_indices()
#        buf[indices] = data.chromosome.ravel()
#        # lines[:, :chromosome_lens] = data.chromosome.ravel()
#        buf[lines.shape.starts+chromosome_lens] = ord("\t")
#        # lines[:, chromosome_lens] = ord("\t")
#
#        buf[line_ends-(end_lens+2)] = ord("\t")
#        buf[line_ends-1] = ord("\n")
#        return buf
#
#    get_data = get_intervals
#
#
#class _VCFBuffer(DelimitedBuffer):
#    dataclass = Variant
#
#    def get_variants(self, fixed_length=False) -> Variant:
#        """Extract variants from VCFBuffer
#
#        Fetches the basic data for a variant from a VCFBuffer and returns a Variant dataset
#
#        Parameters
#        ----------
#        fixed_length : False
#            Wheter or not all sequences are the same length
#
#        Returns
#        -------
#        Variant
#            Variant dataset
#
#        Examples
#        --------
#        5
#
#        """
#        self.validate_if_not()
#        chromosomes = self.get_text(0, fixed_length=False)
#        position = self.get_integers(1).ravel() - 1
#        from_seq = self.get_text(3, fixed_length=fixed_length)
#        to_seq = self.get_text(4, fixed_length=fixed_length)
#        return Variant(chromosomes, position, from_seq, to_seq)
#
#    def get_snps(self):
#        return self.get_variants(fixed_length=True)
#
#    get_data = get_variants
#
#    @classmethod
#    def from_data(cls, data):
#        position = data.position+1
#        position_lens = np.log10(position).astype(int)+1
#        chromosome_lens = data.chromosome.shape.lengths
#
#        ref_lens = 1 # data.ref_seq.shape[-1]
#        alt_lens = 1 # data.alt_seq.shape[-1]
#        line_lengths = chromosome_lens + 1 + position_lens + 1 + 2 + ref_lens + 1 + alt_lens + 1
#        line_ends = np.cumsum(line_lengths)
#        buf = np.empty(line_ends[-1], dtype=np.uint8)
#        lines = RaggedArray(buf, line_lengths)
#        obj = cls(buf, line_ends-1)
#        obj._move_2d_array_to_intervals(cls._move_ints_to_digit_array(position, np.max(position_lens)),
#                                        line_ends-3-ref_lens-alt_lens-2-position_lens, line_ends-3-ref_lens-alt_lens-2)
#
#        indices, _ = RaggedView(lines.shape.starts, chromosome_lens).get_flat_indices()
#        buf[indices] = data.chromosome.ravel()
#
#        ref_indices = lines.shape.starts+chromosome_lens+1+position_lens+2+1
#        alt_indices = lines.shape.starts+chromosome_lens+1+position_lens+2+1+ref_lens+1
#        buf[ref_indices] = data.ref_seq.ravel()
#        buf[alt_indices] = data.alt_seq.ravel()
#        buf[lines.shape.starts+chromosome_lens] = "\t"
#        buf[lines.shape.starts+chromosome_lens+1+position_lens] = "\t"
#        buf[lines.shape.starts+chromosome_lens+1+position_lens+1] = "."
#        buf[lines.shape.starts+chromosome_lens+1+position_lens+2] = "\t"
#        buf[lines.shape.starts+chromosome_lens+1+position_lens+2+1+ref_lens] = "\t"
#        # buf[line_ends-(end_lens+2)] = ord("\t")
#        buf[line_ends-1] = "\n"
#        return buf
'''
