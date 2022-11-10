from bionumpy.io.motifs import read_motif
from bionumpy.simulate.chipseq import simulate_chip_seq_reads, simulate_sequence, ChipSeqSimulationSettings
import matplotlib.pyplot as plt
import bionumpy as bnp
from bionumpy.intervals import get_pileup

chromosome_sizes = {"chr1": 1000,
                    "chr2": 2000,
                    "chr3": 3000}


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
    pileup = get_pileup(reads, length)
    plt.plot(pileup.to_array())
    plt.show()
    return reads


def test_simulate(do_plot=False):
    motif = read_motif("example_data/MA0080.1.jaspar")
    settings = ChipSeqSimulationSettings(motif)
    sequences = {name: simulate_sequence("acgt", size) for name, size in chromosome_sizes.items()}
    multistream = bnp.MultiStream(chromosome_sizes, sequences=sequences)
    reads = simulate_chip_seq_reads(multistream.sequences, settings)
    if do_plot:
        reads = plot_pileup(reads, multistream.lengths)
    print(reads)
    bnp.open("example_data/simulated_chip_seq.bed", "w").write(reads)


if __name__ == "__main__":
    test_simulate(True)
