from npstructures import VarLenArray
from .file_buffers import FileBuffer, NEWLINE
from .datatypes import Interval, Variant, VariantWithGenotypes, SequenceEntry
from .sequences import Sequence
from .encodings import DigitEncoding, GenotypeEncoding, BaseEncoding
import dataclasses
import numpy as np


class DelimitedBuffer(FileBuffer):
    DELIMITER = ord("\t")
    COMMENT = ord("#")

    def __init__(self, data, new_lines):
        super().__init__(data, new_lines)
        self._delimiters = np.concatenate(
            ([-1], np.flatnonzero(self._data == self.DELIMITER), self._new_lines)
        )
        self._delimiters.sort(kind="mergesort")

    @classmethod
    def from_raw_buffer(cls, chunk):
        new_lines = np.flatnonzero(chunk == NEWLINE)
        return cls(chunk[: new_lines[-1] + 1], new_lines)

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
        # delimiters = self._delimiters.reshape(-1, self._n_cols)
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
        chromosomes = VarLenArray(Sequence.from_array(self.get_text(0)))
        positions = self.get_integers(cols=[1, 2])
        return Interval(chromosomes, positions[..., 0], positions[..., 1])

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
        chromosomes = VarLenArray(Sequence.from_array(self.get_text(0)))
        position = self.get_integers(1).ravel() - 1
        from_seq = self.get_text(3, fixed_length=fixed_length)
        to_seq = self.get_text(4, fixed_length=fixed_length)
        return Variant(chromosomes, position, from_seq, to_seq)

    def get_snps(self):
        return self.get_variants(fixed_length=True)

    get_data = get_variants


class VCFMatrixBuffer(VCFBuffer):
    dataclass = VariantWithGenotypes

    def get_entries(self, fixed_length=False):
        self.validate_if_not()
        variants = self.get_variants(fixed_length)
        genotypes = self.get_text_range(np.arange(9, self._n_cols), end=3)
        n_samples = self._n_cols - 9
        genotypes = GenotypeEncoding.encode(genotypes.reshape(-1, n_samples, 3))
        return VariantWithGenotypes(
            variants.chromosome,
            variants.position,
            variants.ref_seq,
            variants.alt_seq,
            genotypes,
        )

    get_data = get_entries


class GfaSequenceBuffer(DelimitedBuffer):
    dataclass = SequenceEntry

    def get_sequences(self):
        ids = self.get_text(1, fixed_length=False)
        sequences = self.get_text(col=2, fixed_length=False)
        return SequenceEntry(ids, sequences)

    get_data = get_sequences


def get_bufferclass_for_datatype(_dataclass):
    class DatatypeBuffer(DelimitedBuffer):
        dataclass = _dataclass

        def get_data(self):
            self.validate_if_not()
            columns = []
            for col_number, field in enumerate(dataclasses.fields(self.dataclass)):
                if field.type is None:
                    col = None
                elif field.type == str:
                    col = self.get_text(col_number, fixed_length=False)
                elif field.type == int:
                    col = self.get_integers(col_number)
                else:
                    assert False, field
                columns.append(col)
            n_entries = len(next(col for col in columns if col is not None))
            columns = [c if c is not None else np.empty((n_entries, 0))
                       for c in columns]
            return self.dataclass(*columns)
    DatatypeBuffer.__name__ = _dataclass.__name__+"Buffer"
    DatatypeBuffer.__qualname__ = _dataclass.__qualname__+"Buffer"
    return DatatypeBuffer
