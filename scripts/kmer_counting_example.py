import bionumpy as bnp


def count_kmers(sequence_entries: bnp.EncodedArray) -> bnp.EncodedCounts:
    sequence = bnp.as_encoded_array(sequence_entries, bnp.DNAEncoding)
    kmers = bnp.get_kmers(sequence, k=5)
    return bnp.count_encoded(kmers, axis=None)


def count_all_kmers(input_file: str, output_file: str):
    buffer_type = bnp.TwoLineFastaBuffer if input_file.endswith(('.fa', '.fa.gz')) else None
    stream = bnp.open(input_file, buffer_type=buffer_type).read_chunks()
    kmers = sum(count_kmers(chunk.sequence) for chunk in stream)
    with open(output_file, 'w') as f:
        f.writelines(f'{kmer}\t{count}\n'
                     for kmer, count
                     in sorted(zip(kmers.alphabet, kmers.counts)))


def test():
    count_all_kmers('example_data/big.fq.gz', 'tmp.csv')


if __name__ == "__main__":
    import sys
    count_all_kmers(*sys.argv[1:])