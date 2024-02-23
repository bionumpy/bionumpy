import bionumpy as bnp


def calculate_forbes(chrom_sizes_file: str, filename_a: str, filename_b: str):
    genome = bnp.Genome.from_file(chrom_sizes_file)

    # Read bed files as binary masks over the genome
    a_mask = genome.read_intervals(filename_a).get_mask()
    b_mask = genome.read_intervals(filename_b).get_mask()
    observed_intersection = (b_mask & a_mask).sum()
    expected_intersection = (a_mask.sum() * b_mask.sum()) / genome.size

    return observed_intersection / expected_intersection


def test():
    similarity = calculate_forbes("example_data/hg38.chrom.sizes",
                                  "example_data/ctcf.bed.gz",
                                  "example_data/znf263.bed.gz")

    print(similarity)


if __name__ == "__main__":
    import sys
    calculate_forbes(*sys.argv[1:])
