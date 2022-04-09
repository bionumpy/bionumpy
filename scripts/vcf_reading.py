import sys
from bionumpy import IndexedFasta, bnp_open, chromosome_map
import numpy as np

@chromosome_map(reduction=sum)
def func(variants, reference):
    return np.bincount(reference[variants.position], minlength=256)

if __name__ == "__main__":
    vcf_filename = sys.argv[1]
    reference_filename = sys.argv[2]
    vcf_file = bnp_open(vcf_filename)
    reference = bnp_open(reference_filename, remove_chr=True)
    print(func(vcf_file, reference))

    #for chromosome, variants in vcf_file:
    #     print(func(variants, reference[chromosome]))
        # 
        # ref_seq = reference[chromosome]
        # s = ref_seq[variants.position]
        # print(chromosome, len(variants))
        # print(s[:10])
        # print(variants[:10])
