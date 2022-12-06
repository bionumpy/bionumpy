import numpy as np
from ..encoded_array import EncodedArray, as_encoded_array, EncodedRaggedArray
from .multiline_buffer import FastaIdxBuffer, FastaIdx
from ..datatypes import Interval
from .files import bnp_open
from ..encodings import BaseEncoding


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

    def __init__(self, filename: str):
        self._filename = filename
        self._index = read_index(filename+".fai")
        self._f_obj = open(filename, "rb")

    def get_contig_lengths(self) -> dict:
        """Return a dict of chromosome names to seqeunce lengths

        Returns
        -------
        dict
            chromosome name to sequence length mapping
        """
        return {name: values["lenc"] for name, values in self._index.items()}

    def keys(self):
        return self._index.keys()

    def values(self):
        return (self[key] for key in self.keys())

    def items(self):
        return ((key, self[key]) for key in self.keys())

    def __repr__(self):
        return f"Indexed Fasta File with chromosome sizes: {self.get_contig_lengths()}"

    def __getitem__(self, chromosome: str) -> EncodedArray:
        """Return entire sequence of the given chromosome

        Parameters
        ----------
        chromosome : str
            chromsome name

        Returns
        -------
        EncodedArray
            The sequence for that chromoeme
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

    def get_interval_sequences(self, intervals: Interval) -> EncodedRaggedArray:
        """Get the seqeunces for a set of genomic intervals

        Parameters
        ----------
        intervals : Interval
            Intervals

        Returns
        -------
        EncodedRaggedArray
            Sequences
        """
        sequences = []
        lengths = []
        delete_indices = []
        cur_offset = 0
        for interval in intervals:
            chromosome = interval.chromosome.to_string()
            idx = self._index[chromosome]
            lenb, rlen, lenc = (idx["lenb"], idx["rlen"], idx["lenc"])
            assert interval.stop <= rlen
            start_row = interval.start//lenc
            start_mod = interval.start % lenc
            start_offset = start_row*lenb+start_mod
            stop_row = interval.stop // lenc
            stop_offset = stop_row*lenb+interval.stop % lenc
            self._f_obj.seek(idx["offset"] + start_offset)
            lengths.append(stop_offset-start_offset-(stop_row-start_row))
            sequences.extend(self._f_obj.read(stop_offset-start_offset))
            delete_indices.extend(cur_offset + lenb*(j+1)-1-start_mod for j in range(stop_row-start_row))
            cur_offset += stop_offset-start_offset
        s = np.delete(np.array(sequences, dtype=np.uint8), delete_indices)
        a = EncodedArray(s, BaseEncoding)
        return EncodedRaggedArray(a, lengths)
