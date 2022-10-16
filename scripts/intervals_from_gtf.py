import bionumpy as bnp
from bionumpy.intervals import merge_intervals, sort_intervals
from bionumpy.groupby import groupby


def main(gtf_filename, out_filename):
    entries = bnp.open(gtf_filename).read_chunks()
    grouped = groupby(entries, "chromosome")
    merged = merge_intervals(sort_intervals(grouped))
    with bnp.open(out_filename, "w") as f:
        f.write(merged)


def test():
    main("example_data/small.gtf", "example_data/small_merged.bed")


if __name__ == "__main__":
    test()
