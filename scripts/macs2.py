import dataclasses
from scipy.special import pdtrc
import numpy as np
import logging

from bionumpy.io.bam import BamIntervalBuffer
from bionumpy.arithmetics import get_pileup, get_boolean_mask, merge_intervals
from bionumpy.datatypes import Interval


# import matplotlib.pyplot as plt
import bionumpy as bnp

logging.basicConfig(level=logging.INFO)


@dataclasses.dataclass
class Params:
    fragment_length: int = 150
    p_value_cufoff: float = 0.001


def extend_to_size(intervals, fragment_length, size):
    is_forward = intervals.strand == "+"
    start = np.where(is_forward,
                     intervals.start,
                     np.maximum(intervals.stop-fragment_length, 0))
    # print(intervals, fragment_length, size)
    stop = np.where(is_forward,
                    intervals.stop,
                    np.minimum(intervals.start+fragment_length, size))

    return dataclasses.replace(
        intervals,
        start=start,
        stop=stop)


def get_fragment_pileup(reads, fragment_length, size):
    logging.info("Getting fragment pileup")
    fragments = extend_to_size(reads, fragment_length, size)
    return get_pileup(fragments, size)


def get_control_pileup(reads, size, window_sizes, read_rate):
    mid_points = (reads.start+reads.stop)//2
    pileup = float(read_rate)
    for window_size in window_sizes:
        start = np.maximum(mid_points-window_size//2, 0)
        stop = np.minimum(mid_points+window_size//2, size)
        windows = dataclasses.replace(reads, start=start, stop=stop)
        new_pileup = get_pileup(windows, size)
        pileup = np.maximum(pileup, new_pileup/window_size)
    return pileup


def logsf(count, mu):
    return np.log(pdtrc(count, mu))


def get_p_values(intervals, chrom_size, fragment_length, read_rate):
    intervals.strand = intervals.strand.ravel()
    fragment_pileup = get_fragment_pileup(intervals, fragment_length, chrom_size)
    print(fragment_pileup.to_bedgraph('.'))
    # print(bnp.bedgraph.from_runlength_array(str(intervals.chromosome[0]), fragment_pileup))
    control = fragment_length*get_control_pileup(intervals, chrom_size, [1000, 10000], read_rate)
    p_values = logsf(fragment_pileup, control)
    return p_values


def call_peaks(p_values, p_value_cufoff, min_length, max_gap=30):
    peaks = p_values < np.log(p_value_cufoff)
    peaks = Interval(['.']*len(peaks.starts), peaks.starts, peaks.ends)[peaks.values]
    peaks = merge_intervals(peaks, distance=max_gap)
    peaks = peaks[(peaks.stop-peaks.start) >= min_length]
    return peaks


@bnp.streamable()
def macs2(intervals, chrom_size, fragment_length, read_rate, p_value_cutoff, min_length, max_gap=30):
    if not len(intervals):
        return Interval.empty()
    p_values = get_p_values(intervals, chrom_size,
                            fragment_length, read_rate)
    peaks = call_peaks(p_values, p_value_cutoff, fragment_length)
    return Interval(intervals.chromosome[:len(peaks)], peaks.start, peaks.stop)


def main(filename: str, genome_file: str, fragment_length: int = 150, p_value_cutoff: float = 0.001, outfilename: str = None):
    genome = bnp.open(genome_file, buffer_type=bnp.io.files.ChromosomeSizeBuffer).read()
    genome_size = genome.size.sum()
    chrom_sizes = {str(name): size for name, size in zip(genome.name, genome.size)}
    intervals = bnp.open(filename, buffer_type=bnp.io.delimited_buffers.Bed6Buffer).read_chunks()
    multistream = bnp.MultiStream(chrom_sizes, intervals=intervals)
    n_reads = bnp.count_entries(filename)
    read_rate = n_reads/genome_size
    result = macs2(multistream.intervals, multistream.lengths, fragment_length,
                   read_rate, p_value_cutoff, fragment_length)
    if outfilename is not None:
        with bnp.open(outfilename, 'w') as f:
            f.write(result)
    return result


def test():
    #res = main("example_data/small_interval.bed", "example_data/small_genome.fa.fai")
    res = main('example_data/simulated_chip_seq.bed', 'example_data/simulated.chrom.sizes', fragment_length=200)
    for chunk in res:
        print(chunk)


def big():
    with bnp.open("/home/knut/Data/peaks.bed", "w") as outfile:
        for data in main("/home/knut/Data/subset.bed", "/home/knut/Data/hg38.chrom.sizes"):
            if len(data):
                outfile.write(data)


if __name__ == "__main__":
    import typer
    typer.run(main)
    #big()
#intervals = bnp.open("/home/knut/Data/ENCFF296OGN.bed", buffer_type=Bed6Buffer).read_chunks()


# intervals = bnp.open("/home/knut/Data/ENCFF296OGN.bam", buffer_type=BamIntervalBuffer).read_chunks(chunk_size=524288*32)
# # intervals = alignment_to_interval(reads)
# # first_chunk = next(chunk_lines(iter(intervals), 10000000))
# grouped = groupby(intervals, "chromosome")
# chrom, first_chunk = next(iter(grouped))
# chrom_size = 248956422
# get_p_values(first_chunk, 300, chrom_size)
