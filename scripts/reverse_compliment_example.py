import bionumpy as bnp


def reverse_complement(input_filename: str, output_filename: str):
    """Reverse complements a fasta or fastq file and writes the result to a new file."""
    bt = lambda filename: (bnp.TwoLineFastaBuffer if filename.endswith(('fa', 'fa.gz')) else None)
    with bnp.open(output_filename, "w", buffer_type=bt(output_filename)) as outfile:
        for chunk in bnp.open(input_filename, buffer_type=bt(input_filename)).read_chunks():
            rc = bnp.sequence.get_reverse_complement(chunk.sequence)
            outfile.write(bnp.replace(chunk, sequence=rc))


def test():
    reverse_complement('example_data/big.fq.gz', 'example_data/big_rc.fq.gz')
    assert bnp.count_entries('example_data/big_rc.fq.gz') == bnp.count_entries('example_data/big.fq.gz')


if __name__ == '__main__':
    import sys
    reverse_complement(sys.argv[1], sys.argv[2])
