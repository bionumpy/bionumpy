import numpy as np

import bionumpy as bnp


def method_name(infile, outfile):
    with bnp.open(outfile, 'w') as f:
        for chunk in bnp.open(infile).read_chunks():
            mask = (chunk.info.LoF.lengths > 0) & chunk.info.ONCOGENE
            print(np.sum(chunk.info.LoF.lengths > 0), chunk.info.ONCOGENE.sum(),
                  mask.sum())

            # print(mask.sum())
            subset = chunk[mask]
            print(subset)
            f.write(subset)
    print('>', bnp.count_entries(outfile))


def test():
    method_name('/home/knut/Downloads/CPCT02020719T_pcgr (1).vcf', 'tmp.vcf')


if __name__ == "__main__":
    import sys

    method_name(sys.argv[1], sys.argv[2])
