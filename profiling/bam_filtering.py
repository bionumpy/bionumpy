import bionumpy as bnp


def filter_bam(input_file, output_file):
    with bnp.open(output_file, 'w') as f:
        for chunk in bnp.open(input_file).read_chunks(min_chunk_size=100000):
            f.write(chunk[chunk.mapq == 60])

def _test_profile():
    filter_bam('../benchmarks/results/bams/big.bam', 'tmp.bam')

