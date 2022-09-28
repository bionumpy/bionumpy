import itertools
import logging
from npstructures import RaggedArray, RaggedView
from .file_buffers import FileBuffer, NEWLINE
from .datatypes import Interval, Variant, VariantWithGenotypes, SequenceEntry
from .sequences import Sequence
import dataclasses
from .encodings import DigitEncoding, GenotypeEncoding, PhasedGenotypeEncoding
import numpy as np


class DelimitedBuffer(FileBuffer):
    """Base class for file buffers for delimited files such as csv or tsv.
    Each line should correspond to an entry, and each column to a variable.

    Provides convenience methods for extracting and decoding integers from columns,
    and text from columns into Sequences objects
    """

    DELIMITER = ord("\t")
    COMMENT = ord("#")

    def __init__(self, data, new_lines, delimiters=None, header_data=None):
        super().__init__(data, new_lines)
        if delimiters is None:
            delimiters = np.concatenate(
                ([-1], np.flatnonzero(self._data == self.DELIMITER), self._new_lines)
            )
            delimiters.sort(kind="mergesort")
        self._delimiters = delimiters
        self._header_data = header_data

    @classmethod
    def from_raw_buffer(cls, chunk, header_data=None):
        mask = chunk == NEWLINE
        mask |= chunk == cls.DELIMITER
        delimiters = np.flatnonzero(mask)
        n_fields = next((i+1 for i, v in enumerate(delimiters) if chunk[v] == NEWLINE), None)
        if n_fields is None:
            logging.warning("Foud no new lines. Chunk size may be too low. Try increasing")
            raise
        new_lines = delimiters[(n_fields-1)::n_fields]
        delimiters = np.concatenate(([-1], delimiters[:n_fields*len(new_lines)]))
        return cls(chunk[:new_lines[-1] + 1], new_lines, delimiters, header_data)

    def get_integers(self, cols) -> np.ndarray:
        """Get integers from integer string

        Extract integers from the specified columns

        Parameters
        ----------
        cols : list
            list of columns containing integers

        Examples
        --------
        FIXME: Add docs.

        """
        cols = np.asanyarray(cols)
        integer_starts = self._delimiters[:-1].reshape(-1, self._n_cols)[:, cols] + 1
        integer_ends = self._delimiters[1:].reshape(-1, self._n_cols)[:, cols]
        integers = self._extract_integers(integer_starts.ravel(), integer_ends.ravel())
        return integers.reshape(-1, cols.size)

    def get_text(self, col, fixed_length=True):
        """Extract text from a column

        Extract strings from the specified column into either a 2d
        array or a RaggedArray

        Parameters
        ----------
        col : int
            column index
        fixed_length : bool
            whether all strings have equal length

        Examples
        --------
        FIXME: Add docs.

        """
        self.validate_if_not()
        starts = self._delimiters[:-1].reshape(-1, self._n_cols)[:, col] + 1
        ends = self._delimiters[1:].reshape(-1, self._n_cols)[:, col]
        if fixed_length:
            return self._move_intervals_to_2d_array(starts, ends)
        else:
            return self._move_intervals_to_ragged_array(starts, ends)

    def get_text_range(self, col, start=0, end=None) -> np.ndarray:
        """Get substrings of a column

        Extract the text from start to end of each entry in column

        Parameters
        ----------
        col : int
            column index
        start : int
            start of substring
        end : int
            end of substring

        Returns
        -------
        np.ndarray
            array containing the extracted substrings

        Examples
        --------
        FIXME: Add docs.

        """
        self.validate_if_not()
        starts = self._delimiters[:-1].reshape(-1, self._n_cols)[:, col] + 1 + start
        if end is not None:
            return self._data[starts[..., np.newaxis] + np.arange(end-start)].reshape(-1, end-start)
        else:
            ends = self._delimiters[1:].reshape(-1, self._n_cols)[:, col]
        return self._move_intervals_to_2d_array(starts.ravel(), ends.ravel())

    def _extract_integers(self, integer_starts, integer_ends):
        digit_chars = self._move_intervals_to_2d_array(
            integer_starts, integer_ends, DigitEncoding.MIN_CODE
        )
        n_digits = digit_chars.shape[-1]
        powers = np.uint32(10) ** np.arange(n_digits)[::-1]
        return DigitEncoding.encode(digit_chars) @ powers

    @staticmethod
    def _move_ints_to_digit_array(ints, n_digits):
        powers = np.uint8(10)**np.arange(n_digits)[::-1]
        ret = (ints[..., None]//powers) % 10
        return ret + ord("0")

    def _validate(self):
        chunk = self._data
        delimiters = self._delimiters[1:]
        n_delimiters_per_line = (
            next(i for i, d in enumerate(delimiters) if chunk[d] == NEWLINE) + 1
        )
        self._n_cols = n_delimiters_per_line
        last_new_line = next(
            i for i, d in enumerate(delimiters[::-1]) if chunk[d] == NEWLINE
        )
        delimiters = delimiters[: delimiters.size - last_new_line]
        assert (
            delimiters.size % n_delimiters_per_line == 0
        ), f"irregular number of delimiters per line ({delimiters.size}, {n_delimiters_per_line})"
        delimiters = delimiters.reshape(-1, n_delimiters_per_line)
        assert np.all(chunk[delimiters[:, -1]] == NEWLINE)
        self._validated = True


class BedBuffer(DelimitedBuffer):
    dataclass = Interval

    def get_intervals(self):
        self.validate_if_not()
        chromosomes = self.get_text(0)
        positions = self.get_integers(cols=[1, 2])
        return Interval(chromosomes, positions[..., 0], positions[..., 1])


    @classmethod
    def from_data(cls, data):
        start_lens = np.log10(data.start).astype(int)+1
        end_lens = np.log10(data.end).astype(int)+1
        chromosome_lens = data.chromosome.shape.lengths
        line_lengths = chromosome_lens + 1 + start_lens + 1 + end_lens + 1
        line_ends = np.cumsum(line_lengths)
        buf = np.empty(line_ends[-1], dtype=np.uint8)
        lines = RaggedArray(buf, line_lengths)
        obj = cls(buf, line_ends-1)
        obj._move_2d_array_to_intervals(cls._move_ints_to_digit_array(data.end, np.max(end_lens)),
                                        line_ends-1-end_lens, line_ends-1)

        obj._move_2d_array_to_intervals(cls._move_ints_to_digit_array(data.start, np.max(start_lens)),
                                        line_ends-2-end_lens-start_lens, line_ends-2-end_lens)

        indices, _ = RaggedView(lines.shape.starts, chromosome_lens).get_flat_indices()
        buf[indices] = data.chromosome.ravel()
        # lines[:, :chromosome_lens] = data.chromosome.ravel()
        buf[lines.shape.starts+chromosome_lens] = ord("\t")
        # lines[:, chromosome_lens] = ord("\t")

        buf[line_ends-(end_lens+2)] = ord("\t")
        buf[line_ends-1] = ord("\n")
        return buf

    get_data = get_intervals


class VCFBuffer(DelimitedBuffer):
    dataclass = Variant

    def get_variants(self, fixed_length=False) -> Variant:
        """Extract variants from VCFBuffer

        Fetches the basic data for a variant from a VCFBuffer and returns a Variant dataset

        Parameters
        ----------
        fixed_length : False
            Wheter or not all sequences are the same length

        Returns
        -------
        Variant
            Variant dataset

        Examples
        --------
        5

        """
        self.validate_if_not()
        chromosomes = self.get_text(0, fixed_length=False)
        position = self.get_integers(1).ravel() - 1
        from_seq = self.get_text(3, fixed_length=fixed_length)
        to_seq = self.get_text(4, fixed_length=fixed_length)
        return Variant(chromosomes, position, from_seq, to_seq)

    def get_snps(self):
        return self.get_variants(fixed_length=True)

    get_data = get_variants

    @classmethod
    def from_data(cls, data):
        position = data.position+1
        position_lens = np.log10(position).astype(int)+1
        chromosome_lens = data.chromosome.shape.lengths

        ref_lens = 1 # data.ref_seq.shape[-1]
        alt_lens = 1 # data.alt_seq.shape[-1]
        line_lengths = chromosome_lens + 1 + position_lens + 1 + 2 + ref_lens + 1 + alt_lens + 1
        line_ends = np.cumsum(line_lengths)
        buf = np.empty(line_ends[-1], dtype=np.uint8)
        lines = RaggedArray(buf, line_lengths)
        obj = cls(buf, line_ends-1)
        obj._move_2d_array_to_intervals(cls._move_ints_to_digit_array(position, np.max(position_lens)),
                                        line_ends-3-ref_lens-alt_lens-2-position_lens, line_ends-3-ref_lens-alt_lens-2)

        indices, _ = RaggedView(lines.shape.starts, chromosome_lens).get_flat_indices()
        buf[indices] = data.chromosome.ravel()

        ref_indices = lines.shape.starts+chromosome_lens+1+position_lens+2+1
        alt_indices = lines.shape.starts+chromosome_lens+1+position_lens+2+1+ref_lens+1
        buf[ref_indices] = data.ref_seq.ravel()
        buf[alt_indices] = data.alt_seq.ravel()
        buf[lines.shape.starts+chromosome_lens] = ord("\t")
        buf[lines.shape.starts+chromosome_lens+1+position_lens] = ord("\t")
        buf[lines.shape.starts+chromosome_lens+1+position_lens+1] = ord(".")
        buf[lines.shape.starts+chromosome_lens+1+position_lens+2] = ord("\t")
        buf[lines.shape.starts+chromosome_lens+1+position_lens+2+1+ref_lens] = ord("\t")
        # buf[line_ends-(end_lens+2)] = ord("\t")
        buf[line_ends-1] = ord("\n")
        return buf


class VCFMatrixBuffer(VCFBuffer):
    dataclass = VariantWithGenotypes
    genotype_encoding = GenotypeEncoding

    def get_entries(self, fixed_length=False):
        self.validate_if_not()
        variants = self.get_variants(fixed_length)
        genotypes = self.get_text_range(np.arange(9, self._n_cols), end=3)
        n_samples = self._n_cols - 9
        genotypes = self.genotype_encoding.from_bytes(genotypes.reshape(-1, n_samples, 3))
        return VariantWithGenotypes(
            variants.chromosome,
            variants.position,
            variants.ref_seq,
            variants.alt_seq,
            genotypes,
        )

    get_data = get_entries


class PhasedVCFMatrixBuffer(VCFMatrixBuffer):
    genotype_encoding = PhasedGenotypeEncoding


class GfaSequenceBuffer(DelimitedBuffer):
    dataclass = SequenceEntry

    def get_sequences(self):
        ids = self.get_text(1, fixed_length=False)
        sequences = self.get_text(col=2, fixed_length=False)
        return SequenceEntry(ids, sequences)

    get_data = get_sequences


def get_bufferclass_for_datatype(_dataclass, delimiter="\t", has_header=False, comment="#"):
    """
    Create a FileBuffer subclass that handles reding of text and integers from a
    buffer in to an @npdataclas object
    """

    class DatatypeBuffer(DelimitedBuffer):
        DELIMITER = ord(delimiter)
        COMMENT = ord(comment)
        dataclass = _dataclass
        fields = None

        def __init__(self, data, new_lines, delimiters=None, header_data=None):
            super().__init__(data, new_lines, delimiters, header_data)
            self.set_fields_from_header(header_data)

        @classmethod
        def read_header(cls, file_object):
            if not has_header:
                return None
            return file_object.readline().decode('ascii').strip().split(chr(cls.DELIMITER))

        def set_fields_from_header(self, columns):
            if not has_header:
                return None
            fields = dataclasses.fields(self.dataclass)
            self.fields = [next(field for field in fields if field.name == col) for col in columns]
            assert np.array_equal(columns, [field.name for field in self.fields])

        def get_data(self):
            self.validate_if_not()
            columns = {}
            fields = self.fields if self.fields is not None else dataclasses.fields(self.dataclass)
            for col_number, field in enumerate(fields):
                if field.type is None:
                    col = None
                elif field.type == str:
                    col = self.get_text(col_number, fixed_length=False)
                elif field.type == int:
                    col = self.get_integers(col_number)
                else:
                    assert False, field
                columns[field.name] = col
            n_entries = len(next(col for col in columns if col is not None))
            columns = {c: value if c is not None else np.empty((n_entries, 0))
                       for c, value in columns.items()}
            return self.dataclass(**columns)
    DatatypeBuffer.__name__ = _dataclass.__name__+"Buffer"
    DatatypeBuffer.__qualname__ = _dataclass.__qualname__+"Buffer"
    return DatatypeBuffer
