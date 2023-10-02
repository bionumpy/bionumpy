from bionumpy.io.motifs import read_motif
from bionumpy.simulate.chipseq import simulate_chip_seq_reads, simulate_sequence, ChipSeqSimulationSettings
import bionumpy as bnp
from bionumpy.arithmetics import get_pileup
from bionumpy.datatypes import ChromosomeSize

chromosome_sizes = {"chr1": 1000,
                    "chr2": 2000,
                    "chr3": 3000}

chrom_sizes = ChromosomeSize(list(chromosome_sizes.keys()),
                             list(chromosome_sizes.values()))


def pass_through(func):
    def call_and_return(func, obj, *args, **kwargs):
        func(obj, *args, **kwargs)
        return obj

    def new_func(stream, *args, **kwargs):
        cls = stream.__class__
        return cls(call_and_return(func, obj, *args, **kwargs)
                   for obj in stream)


@bnp.streamable()
def plot_pileup(reads, length):
    import matplotlib.pyplot as plt
    pileup = get_pileup(reads, length)
    plt.plot(pileup.to_array())
    plt.show()
    return reads


def test_simulate(do_plot=False):
    motif = read_motif("example_data/MA0080.1.jaspar")
    settings = ChipSeqSimulationSettings(motif)
    sequences = {name: simulate_sequence("acgt", size) for name, size in chromosome_sizes.items()}
    multistream = bnp.MultiStream(chromosome_sizes, sequences=sequences)
    reads = simulate_chip_seq_reads(multistream.sequences, settings, multistream.sequence_names)
    if do_plot:
        reads = plot_pileup(reads, multistream.lengths)
    print(reads)
    with bnp.open("example_data/simulated_chip_seq.bed", "w") as f:
        f.write(reads)
    with bnp.open("example_data/simulated.chrom.sizes", "w") as f:
        f.write(chrom_sizes)


if __name__ == "__main__":
    test_simulate(True)
