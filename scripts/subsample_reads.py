import numpy as np
import bionumpy as bnp
snakemake = None


def subsample_reads(input_filename, output_filename):
    n_entries = bnp.count_entries(input_filename, buffer_type=bnp.TwoLineFastaBuffer)
    n_to_subsample = n_entries // 2
    n_entries_left = n_entries
    n_samples_left = n_to_subsample
    out_file = bnp.open(output_filename, "w", buffer_type=bnp.TwoLineFastaBuffer)
    for chunk in bnp.open(input_filename, buffer_type=bnp.TwoLineFastaBuffer).read_chunks():
        if n_entries_left > len(chunk):
            n_samples_in_this_chunk = np.random.hypergeometric(len(chunk), n_entries_left - len(chunk), n_samples_left)
        else:
            n_samples_in_this_chunk = n_samples_left
        rows = np.random.choice(np.arange(len(chunk)), n_samples_in_this_chunk, replace=False)
        n_samples_left -= n_samples_in_this_chunk
        n_entries_left -= len(chunk)
        out_file.write(chunk[rows])
    assert n_samples_left == 0, n_samples_left
    assert n_entries_left == 0, n_entries_left
    out_file.close()

def test_count_entries():
    filename = 'length150_nreads5000000.fa'
    bnp.count_entries(f'benchmarks/results/dna_sequences/{filename}', buffer_type=bnp.TwoLineFastaBuffer)


def test():
    # filename = 'ENCFF689IPX.fq.gz'
    # filename = 'length150_nreads500000.fa'
    filename = 'length150_nreads5000000.fa'
    subsample_reads(f'benchmarks/results/dna_sequences/{filename}',
                    'tmp.fa')
