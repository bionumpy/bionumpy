import bionumpy as bnp
from bionumpy.arithmetics.similarity_measures import forbes, jaccard
from itertools import groupby
import logging
# logging.basicConfig(level="INFO")
from collections import defaultdict


filenames = {"CREM": "ENCFF324ELP.bed.gz",
             "ZNF263": "ENCFF295XBK.bed.gz",
             "FOXA1": "ENCFF497OQD.bed.gz",
             "NR3C1": "ENCFF201BGD.bed.gz"}


base_directory = "/home/knut/Downloads/"

scores = {}
jaccards = {}

chrom_sizes = bnp.open("/home/knut/Data/hg38.chrom.sizes").read()
key_function = {name.to_string(): i for i, name in enumerate(chrom_sizes.name)}.__getitem__
ctcf_peaks = bnp.arithmetics.sort_all_intervals(
    bnp.open("/home/knut/Downloads/ENCFF843VHC.bed.gz").read(),
    chromosome_key_function=key_function)
print(ctcf_peaks)
print([g[0] for g in groupby(ctcf_peaks.chromosome, key=lambda x: x.to_string())])
print(chrom_sizes)

for name, filename in filenames.items():
    tf_peaks = bnp.arithmetics.sort_all_intervals(
        bnp.open(base_directory+filename).read(),
        chromosome_key_function=key_function)

    scores[name] = forbes(chrom_sizes,
                          ctcf_peaks,
                          tf_peaks)
    jaccards[name] = jaccard(chrom_sizes, ctcf_peaks, tf_peaks)

print(scores)
print(jaccards)
