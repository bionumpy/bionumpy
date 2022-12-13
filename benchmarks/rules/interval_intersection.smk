import bionumpy as bnp


rule intersect_bionumpy:
    input:
        a="results/intervals/{b}.bed",
        b="results/intervals/{a}.bed",
        chrom_sizes="results/chrom.sizes"
    output:
        "results/bionumpy/intersect/{a}_{b}.bed"
    benchmark:
        "benchmarks/intersect/bionumpy/{a}_{b}.txt"
    run:
        chrom_sizes = bnp.open(input.chrom_sizes).read()
        a = bnp.open(input[0], buffer_type=bnp.Bed6Buffer).read_chunks()
        b = bnp.open(input[1], buffer_type=bnp.Bed6Buffer).read_chunks()
        ms = bnp.streams.MultiStream(chrom_sizes, a=a, b=b)
        result = bnp.arithmetics.intersect(ms.a, ms.b)
        with bnp.open(output[0], "w", buffer_type=bnp.Bed6Buffer) as f:
            f.write(result)


rule intersect_bedtools:
    input:
        left="results/intervals/{b}.bed",
        right="results/intervals/{a}.bed"
    output:
        "results/bedtools/intersect/{a}_{b}.bed"
    benchmark:
        "benchmarks/intersect/bedtools/{a}_{b}.txt"
    log:
        "logs/bedtools/intersect/{a}_{b}.log"
    params:
        # bedtools intersect is faster when it knows input is sorted
        extra="-sorted"
    wrapper:
        "v1.19.2/bio/bedtools/intersect"
