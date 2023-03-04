import bionumpy as bnp
from bionumpy.genomic_data.global_offset import GlobalOffset
import time
import dataclasses


@dataclasses.dataclass
class Input:
    chrom_sizes: str
    a: str
    b: str

@dataclasses.dataclass
class Snakefile:
    input: Input
    output: list

try:
    snakemake
except:
    snakemake = Snakefile(Input('../../example_data/hg38.chrom.sizes',
                                '../results/intervals/ENCFF227NIG_mapped_reads_1m.bed',
                                '../results/bed_files/ENCFF491EEI.bed.gz'),
                          ['../../tmp.bed'])
                          
t = time.perf_counter()
chrom_sizes = bnp.open(snakemake.input.chrom_sizes).read()

global_offset = GlobalOffset(chrom_sizes)
a = bnp.open(snakemake.input.a, buffer_type=bnp.Bed6Buffer).read()
b = bnp.open(snakemake.input.b, buffer_type=bnp.Bed6Buffer).read()

a = global_offset.from_local_interval(a)
b = global_offset.from_local_interval(b)

#chromosome_encoding = bnp.encodings.string_encodings.StringEncoding(chrom_sizes.name)
#a = dataclasses.replace(a, chromosome=chromosome_encoding.encode(a.chromosome))
# b = dataclasses.replace(b, chromosome=chromosome_encoding.encode(b.chromosome))
result = bnp.arithmetics.unique_intersect(a, b, sum(chrom_sizes.size))
result = global_offset.to_local_interval(result)
with bnp.open(snakemake.output[0], "w", buffer_type=bnp.Bed6Buffer) as f:
    f.write(result)
print(time.perf_counter() - t)
