
rule bionumpy_sequence_length_distribution:
    input:
        "results/dna_sequences/{filename}.fq.gz"
    output:
        "results/bionumpy/sequence_length_distribution/{filename}.csv"
    benchmark:
        "benchmarks/sequence_length_distribution/bionumpy/{filename}.txt"
    shell:
        "python ../scripts/sequence_length_distribution_example.py {input} {output}"
    
rule awk_sequence_length_distribution:
    input:
        "results/dna_sequences/{filename}.fq.gz"
    output:
        "results/awk/sequence_length_distribution/{filename}.csv"
    benchmark:
        "benchmarks/sequence_length_distribution/awk/{filename}.txt"
    conda:
        "../envs/awk.yml"
    shell:
        "zcat {input} | awk 'NR%4 == 2 {{lengths[length($0)]++}} END {{for (l in lengths) {{print l, lengths[l]}}}}' > {output}"         



rule python_sequence_length_distribution:
    input:
        "results/dna_sequences/{filename}.fq.gz"
    output:
        "results/python/sequence_length_distribution/{filename}.csv"
    benchmark:
        "benchmarks/sequence_length_distribution/python/{filename}.txt"
    run:
        import gzip
        from collections import defaultdict
        counts = defaultdict(int)
        with gzip.open(input[0]) as f:
            for i, line in enumerate(f):
                line = line.decode("utf-8")
                if i % 4 != 1:
                    continue

                length = len(line.strip())
                counts[length] += 1

        sorted_counts = sorted(counts)

        with open(output[0], "w") as f:
            f.write("\n".join(
                ["%d %d" % (c, counts[c]) for c in sorted_counts]
            ) + "\n")



