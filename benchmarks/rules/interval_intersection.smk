import bionumpy as bnp


rule intersect_bionumpy:
    input:
        a="results/intervals/{b}.bed",
        b="results/intervals/{a}.bed",
        chrom_sizes="../example_data/hg38.chrom.sizes"
    output:
        "results/bionumpy/intersect/{a}-vs-{b}.bed"
    benchmark:
        "benchmarks/intersect/bionumpy/{a}-vs-{b}.txt"
    run:
        chrom_sizes = bnp.open(input.chrom_sizes).read()
        a = bnp.open(input[0], buffer_type=bnp.Bed6Buffer).read()
        b = bnp.open(input[1], buffer_type=bnp.Bed6Buffer).read()
        print("Sorting")
        a = bnp.arithmetics.sort_intervals(a, sort_order=chrom_sizes.name.tolist())
        b = bnp.arithmetics.sort_intervals(b, sort_order=chrom_sizes.name.tolist())
        print("DOne sorting")

        #ms = bnp.streams.MultiStream(chrom_sizes, a=a, b=b)
        result = bnp.arithmetics.intersect(a, b)
        with bnp.open(output[0], "w", buffer_type=bnp.Bed6Buffer) as f:
            f.write(result)


rule intersect_bedtools:
    input:
        left="results/intervals/{b}.bed",
        right="results/intervals/{a}.bed"
    output:
        "results/bedtools/intersect/{a}-vs-{b}.bed"
    benchmark:
        "benchmarks/intersect/bedtools/{a}-vs-{b}.txt"
    log:
        "logs/bedtools/intersect/{a}-vs-{b}.log"
    params:
        # bedtools intersect is faster when it knows input is sorted
        extra=""
    wrapper:
        "v1.19.2/bio/bedtools/intersect"
