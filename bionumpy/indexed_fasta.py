from .chromosome_provider import ChromosomeDictProvider
import numpy as np

from pyfaidx import Fasta, Faidx


class IndexedFasta(ChromosomeDictProvider):
    def __init__(self, filename, add_chr=False, remove_chr=False):
        self._filename = filename
        self._index = Faidx(filename).index
        self._f_obj = open(filename, "rb")
        self._add_chr = add_chr
        self._remove_chr = remove_chr

    def __getitem__(self, chromosome):
        if self._add_chr:
            assert chromosome.startswith("chr"), chromosome
            chromosome=chromosome[3:]
        if self._remove_chr:
            assert not chromosome.startswith("chr"), chromosome
            chromosome="chr"+chromosome
        if chromosome not in self._index:
            print(chromosome, list(self._index.keys()))
        idx = self._index[chromosome]
        print(idx)
        lenb, rlen, lenc = (idx["lenb"], idx["rlen"], idx["lenc"])
        n_rows = (rlen+lenc-1)//lenc
        data = np.empty(lenb*n_rows, dtype=np.uint8)

        self._f_obj.seek(idx.offset)
        self._f_obj.readinto(data[:rlen])
        data = data.reshape(n_rows, lenb)
        ret = data[:, :idx["lenc"]].ravel()[:idx["rlen"]]
        assert ret.size == idx["rlen"], (ret.size, idx["rlen"], ret.size -idx["rlen"], data.shape)
        return ret
    
to_str = lambda x: "".join(chr(c) for c in x)

if __name__ == "__main__":

    filename = "/home/knut/Data/human_g1k_v37.fasta"
    f = Fasta(filename)
    f2 = IndexedFasta(filename, add_chr=True)
    print(f["2"][100000:100020])
    print(f2["chr2"][100000:100020])
    print(to_str(f2["chr2"][100000:100020]))
