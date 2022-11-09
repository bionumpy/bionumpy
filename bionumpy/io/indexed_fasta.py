from ..encoded_array import EncodedArray
from ..streams.grouped import grouped_dict

import numpy as np


def read_index(filename):
    split_lines = (line.split() for line in open(filename))
    return {chromosome:
            {"rlen": int(rlen), "offset": int(offset), "lenc": int(lenc), "lenb": int(lenb)}
            for chromosome, rlen, offset, lenc, lenb in split_lines}


class IndexedFasta:
    def __init__(self, filename, add_chr=False):
        self._filename = filename
        self._index = read_index(filename+".fai")# Faidx(filename).index
        self._f_obj = open(filename, "rb")

    def get_contig_lengths(self):
        return {name: values["lenc"] for name, values in self._index.items()}

    def __get_chrom_name(self, name):
        if name in self._index:
            return name
        elif name.startswith("chr") and name[3:] in self._index:
            return name[3:]
        else:
            if ("chr"+name) in self._index:
                return "chr" + name
        assert False

    def __getitem__(self, chromosome):
        idx = self._index[self.__get_chrom_name(chromosome)]
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
        return EncodedArray(((ret - ord("A")) % 32) + ord("A"))


to_str = lambda x: "".join(chr(c) for c in x)

if __name__ == "__main__":
    from pyfaidx import Fasta, Faidx

    filename = "/home/knut/Data/human_g1k_v37.fasta"
    f = Fasta(filename)
    f2 = IndexedFasta(filename, add_chr=True)
