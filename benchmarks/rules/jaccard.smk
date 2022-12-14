import itertools


tf_experiments = ["ENCFF843VHC", "ENCFF324ELP", "ENCFF201BGD", "ENCFF295XBK", "ENCFF497OQD", "ENCSR140GLO"]
tf_experiments = [l.strip() for l in open("../example_data/encode_k562_tf_ids.txt")]
sorted(tf_experiments)


def bedtools_commands(wildcards):
    n_tfs = int(wildcards.n_tfs)
    tf_pairs = set(itertools.combinations(tf_experiments[0:n_tfs],2))
    path = "results/bed_files/"
    file_ending = ".sorted.bed"
    # A bit hacky, in order to run all bedtools commands in one rule so benchmark is valid
    bedtools_commands = "\n".join([path + tf1 + file_ending + "," + path + tf2 + file_ending
                               for tf1, tf2 in tf_pairs])
    return bedtools_commands


def bed_input(wildcards):
    n_tfs = int(wildcards.n_tfs)
    return ["results/bed_files/" + tf + ".sorted.bed" for tf in tf_experiments[0:n_tfs]]


rule make_bedtools_input:
    input:
        bed_input
    output:
        "results/bedtools/jaccard_all_vs_all/input{n_tfs}.txt"
    params:
        bedtools_commands=bedtools_commands
    script:
        "../scripts/make_bedtools_input.py"


rule pyranges_jaccard_all_vs_all:
    input:
        bed_input
    output:
        "results/pyranges/jaccard_all_vs_all/ntfs{n_tfs}.txt"
    benchmark:
        "benchmarks/jaccard_all_vs_all/pyranges/ntfs{n_tfs}.txt"
    run:
        import pyranges
        files = {i: pyranges.read_bed(i) for i in input}
        results = {}
        n = 0
        for a in files:
            for b in files:
                run_id = frozenset([a, b])
                if a == b or run_id in results:
                    continue
                j = files[a].stats.jaccard(files[b])
                results[run_id] = j
                n += 1

        print("%d pairs compared" % n)

        with open(output[0], "w") as f:
            f.write(str(list(results.values())) + "\n")



rule bedtools_jaccard_all_vs_all:
    input:
        "results/bedtools/jaccard_all_vs_all/input{n_tfs}.txt"
    output:
        "results/bedtools/jaccard_all_vs_all/ntfs{n_tfs}.txt"
    params:
        bedtools_commands=bedtools_commands
    conda:
        "../envs/bedtools.yml"
    benchmark:
        "benchmarks/jaccard_all_vs_all/bedtools/ntfs{n_tfs}.txt"
    shell:
        """
        touch {output}
        cat {input} | parallel -j 1 --colsep "," bedtools jaccard -a {{1}} -b {{2}} >> {output} 
        """


rule bionumpy_jaccard_all_vs_all:
    input:
        bed_input
    output:
        "results/bionumpy/jaccard_all_vs_all/ntfs{n_tfs}.txt"
    benchmark:
        "benchmarks/jaccard_all_vs_all/bionumpy/ntfs{n_tfs}.txt"
    shell:
        """
        python ../scripts/jaccard_all_vs_all_example.py ../example_data/hg38_unix_sorted.chrom.sizes {input} > {output}
        """


rule bionumpy_jaccard_two_bed_files:
    input:
        a="results/bed_files/{a}.bed.gz",
        b="results/bed_files/{b}.bed.gz"
    output:
        "results/bionumpy/jaccard_two_bed_files/{a}-vs-{b}.txt"
    benchmark:
        "benchmarks/jaccard_two_bed_files/bionumpy/{a}-vs-{b}.txt"
    shell:
        """
        python ../scripts/jaccard_all_vs_all_example.py ../example_data/hg38_unix_sorted.chrom.sizes {input} > {output}
        """

# also does sorting in this rule so that's included in the benchmark
rule bedtools_jaccard_two_bed_files:
    input:
        a="results/bed_files/{a}.bed.gz",
        b="results/bed_files/{b}.bed.gz"
    output:
        "results/bedtools/jaccard_two_bed_files/{a}-vs-{b}.txt"
    benchmark:
        "benchmarks/jaccard_two_bed_files/bedtools/{a}-vs-{b}.txt"
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        zcat {input.a} | sort -k1,1 -k2,2n > {input.a}.sorted.bed.tmp        
        zcat {input.b} | sort -k1,1 -k2,2n > {input.b}.sorted.bed.tmp        
        bedtools jaccard -a {input.a}.sorted.bed.tmp -b {input.b}.sorted.bed.tmp | cut -f 3 | tail -n 1 > {output}
        """



rule pyranges_jaccard_two_bed_files:
    input:
        a="results/bed_files/{a}.bed.gz",
        b="results/bed_files/{b}.bed.gz"
    output:
        "results/pyranges/jaccard_two_bed_files/{a}-vs-{b}.txt"
    benchmark:
        "benchmarks/jaccard_two_bed_files/pyranges/{a}-vs-{b}.txt"
    run:
        import pyranges
        a = pyranges.read_bed(input.a)
        b = pyranges.read_bed(input.b)
        j = a.stats.jaccard(b)
        with open(output[0], "w") as f:
            f.write(str(j) + "\n") 
