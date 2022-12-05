import bionumpy as bnp
from bionumpy.arithmetics.similarity_measures import forbes, jaccard
from itertools import groupby


filenames = {"CREM": "ENCFF324ELP.bed.gz",
             "ZNF263": "ENCFF295XBK.bed.gz",
             "FOXA1": "ENCFF497OQD.bed.gz",
             "NR3C1": "ENCFF201BGD.bed.gz"}


ctcf_filename = "ENCFF843VHC.bed.gz"

base_directory = "/home/knut/Downloads/"

scores = {}
jaccards = {}

chrom_sizes = bnp.open("/home/knut/Data/hg38.chrom.sizes").read()
sort_order = chrom_sizes.name.tolist()

ctcf_peaks = bnp.arithmetics.sort_all_intervals(
    bnp.open(base_directory + ctcf_filename).read(), sort_order=sort_order)

for name, filename in filenames.items():
    tf_peaks = bnp.arithmetics.sort_all_intervals(
        bnp.open(base_directory+filename).read(),
        sort_order=sort_order)

    scores[name] = forbes(chrom_sizes,
                          ctcf_peaks,
                          tf_peaks)
    jaccards[name] = jaccard(chrom_sizes, ctcf_peaks, tf_peaks)

print(scores)
print(jaccards)
