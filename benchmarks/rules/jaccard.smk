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


rule bedtools_jaccard:
    input:
        tf1="results/bed_files/{tf1}.sorted.bed",
        tf2="results/bed_files/{tf2}.sorted.bed"
    output:
        "results/bedtools/jaccard/{tf1}.{tf2}.txt"
    conda:
        "../envs/bedtools.yml"
    shell:
        "bedtools jaccard -a {input.tf1} -b {input.tf2} | cut -f 3 | tail -n 1 > {output}"


rule make_bedtools_input:
    input:
        bed_input
    output:
        "results/bedtools/jaccard_all_vs_all/input{n_tfs}.txt"
    params:
        bedtools_commands=bedtools_commands
    script:
        "../scripts/make_bedtools_input.py"


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
