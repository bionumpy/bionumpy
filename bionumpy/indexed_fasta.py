from .chromosome_map import ChromosomeDictProvider
import numpy as np


class IndexedFasta(ChromosomeDictProvider):
    def __init__(self, filename, add_chr=False):
        self._filename = filename
        self._index = Faidx(filename).index
        self._f_obj = open(filename, "rb")
        self._add_chr = add_chr

    def __getitem__(self, chromosome):
        if self._add_chr:
            assert chromosome.startswith("chr")
            chromosome=chromosome[3:]
        idx = self._index[chromosome]
        lenb = idx["lenb"]
        n_rows = (idx["rlen"]+lenb-1)//lenb
        data = np.empty(lenb*n_rows, dtype=np.uint8)

        self._f_obj.seek(idx.offset)
        self._f_obj.readinto(data[:idx["rlen"]])
        data = data.reshape(n_rows, lenb)
        return data[:, :idx["lenc"]].ravel()[:idx["rlen"]]
      
to_str = lambda x: "".join(chr(c) for c in x)

if __name__ == "__main__":
    from pyfaidx import Fasta, Faidx
    filename = "/home/knut/Data/human_g1k_v37.fasta"
    f = Fasta(filename)
    f2 = IndexedFasta(filename, add_chr=True)
    print(f["2"][100000:100020])
    print(f2["chr2"][100000:100020])
    print(to_str(f2["chr2"][100000:100020]))
