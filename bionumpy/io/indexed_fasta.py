from pathlib import Path
from typing import Union, Dict, Iterable, Tuple

import numpy as np
from ..encoded_array import EncodedArray, as_encoded_array, EncodedRaggedArray
from .multiline_buffer import FastaIdxBuffer, FastaIdx
from ..datatypes import Interval
from .files import bnp_open
from ..encodings import BaseEncoding
from ..encodings.string_encodings import StringEncoding


def read_index(filename: str) -> dict:
    """Read a fa.fai into a nested dict

    Parameters
    ----------
    filename : str
        The filename for the fasta index

    Returns
    -------
    dict
        nested dict. chromosome names to dicts of index values

    """
    split_lines = (line.split("\t") for line in open(filename))
    return {chromosome.split()[0]:
            {"rlen": int(rlen), "offset": int(offset),
             "lenc": int(lenc), "lenb": int(lenb)}
            for chromosome, rlen, offset, lenc, lenb in split_lines}


def create_index(filename: str) -> FastaIdx:
    """Create a fasta index for a fasta file

    Parameters
    ----------
    filename : str
        Filename of the fasta file

    Returns
    -------
    FastaIdx
        Fasta index as bnpdataclass

    """

    reader = bnp_open(filename, buffer_type=FastaIdxBuffer)
    indice_builders = list(reader.read_chunks())
    offsets = np.cumsum([0]+[idx.byte_size[0] for idx in indice_builders])
    return np.concatenate([
        FastaIdx(idx.chromosome,
                 idx.length,
                 idx.start+offset,
                 idx.characters_per_line,
                 idx.line_length)
        for idx, offset in zip(indice_builders, offsets)])


class IndexedFasta:
    """
    Class representing an indexed fasta file.
    Behaves like dict of chrom names to sequences
    """

    def __init__(self, filename: Union[str, Path]):
        if isinstance(filename, str):
            filename = Path(filename)
        self._filename = filename
        self._index = read_index(filename.with_suffix(filename.suffix + ".fai"))
        self._f_obj = open(filename, "rb")
        self._index_table = FastaIdx.from_entry_tuples(
            [(name, var['rlen'], var['offset'], var['lenc'], var['lenb'])
             for name, var in self._index.items()])#  if '_' not in name])

    def get_contig_lengths(self) -> Dict[str, int]:
        """Return a dict of chromosome names to seqeunce lengths

        Returns
        -------
        dict
            chromosome name to sequence length mapping
        """
        return {name: values["lenc"] for name, values in self._index.items()}

    def keys(self) -> Iterable[str]:
        return self._index.keys()

    def values(self) -> Iterable[EncodedArray]:
        return (self[key] for key in self.keys())

    def items(self) -> Iterable[Tuple[str, EncodedArray]]:
        return ((key, self[key]) for key in self.keys())

    def __repr__(self):
        return f"Indexed Fasta File with chromosome sizes: {self.get_contig_lengths()}"

    def __getitem__(self, chromosome: str) -> EncodedArray:
        """Return entire sequence of the given chromosome

        Parameters
        ----------
        chromosome : str
            chromosome name

        Returns
        -------
        EncodedArray
            The sequence for that chromosome
        """
        idx = self._index[chromosome]
        lenb, rlen, lenc = (idx["lenb"], idx["rlen"], idx["lenc"])
        n_rows = (rlen + lenc - 1) // lenc
        data = np.empty(lenb * n_rows, dtype=np.uint8)
        bytes_to_read = (n_rows - 1) * lenb + (rlen - (n_rows - 1) * lenc)
        self._f_obj.seek(idx["offset"])
        self._f_obj.readinto(data[:bytes_to_read])
        assert np.all(data[:bytes_to_read] > 0), data[:bytes_to_read]
        data = data.reshape(n_rows, lenb)
        ret = data[:, :lenc].ravel()[:rlen]
        assert np.all(ret[:rlen] > 0), ret
        assert ret.size == idx["rlen"], (
            ret.size,
            idx["rlen"],
            ret.size - idx["rlen"],
            data.shape,
        )
        return EncodedArray(ret, BaseEncoding)
        return EncodedArray(((ret - ord("A")) % 32) + ord("A"), BaseEncoding)

    def _get_interval_sequences_fast(self, intervals: Interval) -> EncodedRaggedArray:
        tmp = {name: self._index[name] for name in intervals.chromosome.encoding.get_labels()}
        index_table = FastaIdx.from_entry_tuples(
            [(name, var['rlen'], var['offset'], var['lenc'], var['lenb'])
             for name, var in tmp.items()])

        pre_alloc = np.empty((intervals.stop-intervals.start).sum(), dtype=np.uint8)
        chromosome_i = intervals.chromosome.raw()
        indices: FastaIdx = index_table[chromosome_i]
        start_rows = intervals.start//indices.characters_per_line
        start_mods = intervals.start % indices.characters_per_line
        start_offsets = start_rows*indices.line_length+start_mods

        stop_rows = intervals.stop // indices.characters_per_line
        stop_offsets = stop_rows*indices.line_length+intervals.stop % indices.characters_per_line
        read_starts = indices.start + start_offsets
        read_lengths = stop_offsets-start_offsets

        lengths = intervals.stop-intervals.start
        n_rows = stop_rows-start_rows
        offsets = np.insert(np.cumsum(lengths), 0, 0)
        for read_start, read_length, n_row, start_mod, lenb, a_offset in zip(
                read_starts, read_lengths, n_rows, start_mods, indices.line_length, offsets):
            self._f_obj.seek(read_start)
            r_sequence = np.frombuffer(self._f_obj.read(read_length), dtype=np.uint8)
            sequence = np.delete(r_sequence,
                                 [lenb*(j+1)-1-start_mod
                                  for j in range(n_row)])
            pre_alloc[a_offset:a_offset+sequence.size] = sequence
        a = EncodedArray(pre_alloc, BaseEncoding)
        return EncodedRaggedArray(a, lengths)

    def get_interval_sequences(self, intervals: Interval) -> EncodedRaggedArray:
        """Get the sequences for a set of genomic intervals

        Parameters
        ----------
        intervals : Interval
            Intervals

        Returns
        -------
        EncodedRaggedArray
            Sequences
        """
        if isinstance(intervals.chromosome.encoding, StringEncoding):
            return self._get_interval_sequences_fast(intervals)
        lengths = []
        cur_offset = 0
        pre_alloc = np.empty((intervals.stop-intervals.start).sum(), dtype=np.uint8)
        alloc_offset = 0
        
        for interval in intervals:
            chromosome = interval.chromosome.to_string()
            idx = self._index[chromosome]
            lenb, rlen, lenc = (idx["lenb"], idx["rlen"], idx["lenc"])
            start_row = interval.start//lenc
            start_mod = interval.start % lenc
            start_offset = start_row*lenb+start_mod
            stop_row = interval.stop // lenc
            stop_offset = stop_row*lenb+interval.stop % lenc
            self._f_obj.seek(idx["offset"] + start_offset)
            lengths.append(stop_offset-start_offset-(stop_row-start_row))
            D = stop_offset-start_offset
            tmp = np.frombuffer(self._f_obj.read(stop_offset-start_offset),
                                dtype=np.uint8)
            tmp = np.delete(tmp, [lenb*(j+1)-1-start_mod
                                  for j in range(stop_row-start_row)])
            pre_alloc[alloc_offset:alloc_offset+tmp.size] = tmp
            alloc_offset += tmp.size
            cur_offset += stop_offset-start_offset
        assert alloc_offset == pre_alloc.size, (alloc_offset, pre_alloc.size)
        assert np.all(pre_alloc> 0), np.sum(pre_alloc==0)
        a = EncodedArray(pre_alloc, BaseEncoding)
        return EncodedRaggedArray(a, lengths)
