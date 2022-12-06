

rule simulate_dna:
    output:
        "results/dna_sequences/length{length}_nreads{number}.fa"
    run:
        import random
        with open(output[0], "w") as f:
            for i in range(int(wildcards.number)):
                f.write(f">{i}\n")
                f.write("".join(random.choices("ACGT", k=int(wildcards.length)))+ "\n")


rule make_chrom_sizes:
    output:
        "results/chrom.sizes"
    run:
        out = ""
        for chrom in range(1, 11):
            out += f"{chrom}\t1000000\n"
        open(output[0], "w").write(out)


rule simulate_intervals:
    input:
        "results/chrom.sizes"
    output:
        "results/intervals/nintervals{n_intervals}.bed"
    run:
        import random
        chrom_sizes = [(str(l.split()[0]), int(l.split()[1])) for l in open(input[0])]
        with open(output[0], "w") as f:
            for chromosome_name, size in chrom_sizes:
                start = 1
                for _ in range(int(wildcards.n_intervals)):
                    end = start+random.randint(1, 100)
                    f.write(f"{chromosome_name}\t{start}\t{end}\t.\t0\t+\n")
                    start = end + random.randint(1, 10)

