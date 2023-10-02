import numpy as np
import bionumpy as bnp
import plotly.express as px


def div(stream=False):
    # read a genome
    genome = bnp.Genome.from_file("example_data/hg38.chrom.sizes")

    # read peaks as intervals
    peaks = genome.read_intervals("example_data/ctcf_chr21-22.bed.gz", stream=False)

    # Get reads as intervals
    reads = genome.read_intervals("example_data/ctcf_chr21-22.bam", stream=True)

    # reads are just intervals. We can get a pileup across the genome:
    read_pileup = reads.get_pileup()
    print(read_pileup)

    # This pileup can be index, e.g. by our peaks to get a pileup for each peak
    peaks_pileup = read_pileup[peaks]

    # This Pileup can be indexed like a RaggedArray, and we can e.g. get the max pileup
    # value for each peak
    max_pileup_value_per_peak = np.max(peaks_pileup, axis=1)

    # We can filter the peak intervals
    high_peaks = peaks[max_pileup_value_per_peak > 4]
    print(high_peaks.compute())


"""
We have just released a new version that makes it a lot easier to work with genomic intervals 
and data on a reference genome. Here are a few examples to illustrate how BioNumPy can be used.
"""

"""
Example 1: Creating a pileup from a BAM-file. 
"""
def test_example_2_bam_pileup():
    import bionumpy as bnp
    import plotly.express as px

    # Reading a genome and reads from a bam file
    genome = bnp.Genome.from_file("example_data/chr21-22.chrom.sizes")
    reads = genome.read_intervals("example_data/ctcf_chr21-22.bam")

    # Getting read pileup (stored efficiently as a RunLengthArray)
    pileup = reads.get_pileup()

    # We can index any region
    region = pileup["chr22"][19970400:19970800]

    #px.line(region.to_array()).show()



def _example2():
    # Getting read pileup for certain regions
    genome = bnp.Genome.from_file("example_data/hg38.chrom.sizes")
    reads = genome.read_intervals("example_data/ctcf_chr21-22.bam", stream=False).get_pileup()

    ctcf_peaks = genome.read_intervals("example_data/ctcf_chr21-22.bed.gz", stream=False)
    znf_peaks = genome.read_intervals("example_data/myc_chr21-22.bed.gz", stream=False)

    ctcf_pileup = reads[ctcf_peaks]
    znf_pileup = reads[znf_peaks]

    ctcf_scores = np.max(ctcf_pileup, axis=-1)
    znf_scores = np.max(znf_pileup, axis=-1)

    print(np.mean(ctcf_scores), np.mean(znf_scores))
    print(np.std(ctcf_scores), np.std(znf_scores))

    return


def _test_pileup_around_gene_starts():
    genome = bnp.Genome.from_file("example_data/hg38.chrom.sizes")
    genes = genome.read_intervals("example_data/chr21a22.gtf")
    print(genes)


def trash():
    mean_length = int(np.mean(ctcf_peaks.stop-ctcf_peaks.start))
    average_peak_pileup = np.mean(ctcf_pileup.to_array()[mean_length-20:mean_length+20, :], axis=0)
    fig = px.line(average_peak_pileup)
    fig.show()
    return

    import matplotlib.pyplot as plt
    plt.hist([ctcf_scores, znf_scores], bins=30, stacked=True)
    plt.show()
    return

    import plotly.graph_objects as go
    fig = go.Figure()
    fig.add_trace(go.Histogram(x=znf_scores))
    fig.add_trace(go.Histogram(x=ctcf_scores))
    fig.update_layout(barmode='stack')
    fig.show()


def _example3():
    # lower pileup values around variants
    pass


if __name__ == "__main__":
    _example2()
