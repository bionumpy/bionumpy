import dataclasses
from bionumpy.bam import BamIntervalBuffer
from bionumpy.intervals import RawInterval, merge_intervals
from bionumpy.groupby import groupby
from bionumpy.bedgraph import get_pileup
from bionumpy.chromosome_map import ChromosomeMap
from bionumpy.chromosome_provider import GroupedDict
from bionumpy.datatypes import Interval
from scipy.special import pdtrc

# import matplotlib.pyplot as plt
import bionumpy as bnp
import numpy as np
import logging

logging.basicConfig(level=logging.INFO)
# Extend the intervals to be size 150. Starts should


def passthrogh_write(writer, stream):
    for elem in stream:
        writer.write(elem)
        yield elem


@dataclasses.dataclass
class Params:
    fragment_length: int = 150
    p_value_cufoff: float = 0.001


def extend_to_size(intervals, fragment_length, size):
    is_forward = intervals.strand == "+"
    start = np.where(is_forward,
                     intervals.start,
                     np.maximum(intervals.end-fragment_length, 0))
    end = np.where(is_forward,
                   intervals.end,
                   np.minimum(intervals.start+fragment_length, size))

    return dataclasses.replace(
        intervals,
        start=start,
        end=end)


def get_fragment_pileup(reads, fragment_length, size):
    logging.info("Getting fragment pileup")
    fragments = extend_to_size(reads, fragment_length, size)
    fragments = fragments[np.argsort(fragments.start, kind="mergesort")]
    return get_pileup(fragments, size)


def get_control_pileup(reads, size, window_sizes, read_rate):
    mid_points = (reads.start+reads.end)//2
    pileup = read_rate
    for window_size in window_sizes:
        start = np.maximum(mid_points-window_size//2, 0)
        end = np.minimum(mid_points+window_size//2, size)
        windows = dataclasses.replace(reads, start=start, end=end)
        pileup = np.maximum(pileup, get_pileup(windows, size)/window_size)
    return pileup


def logsf(count, mu):
    return np.log(pdtrc(count, mu))


@ChromosomeMap()
def get_p_values(intervals, chrom_size, fragment_length, read_rate):
    intervals.strand = intervals.strand.ravel()
    fragment_pileup = get_fragment_pileup(intervals, fragment_length, chrom_size)
    control = fragment_length*get_control_pileup(intervals, chrom_size, [10000], read_rate)
    p_values = logsf(fragment_pileup, control)
    return p_values


@ChromosomeMap()
def call_peaks(p_values, p_value_cufoff, min_length, max_gap=30):
    peaks = p_values < np.log(p_value_cufoff)
    peaks = RawInterval(peaks.starts, peaks.ends)[peaks.values]
    peaks = merge_intervals(peaks, distance=max_gap)
    peaks = peaks[(peaks.end-peaks.start) >= min_length]
    return peaks


@ChromosomeMap()
def macs2(intervals, chrom_size, fragment_length, read_rate, p_value_cutoff, min_length, max_gap=30):
    p_values = get_p_values(intervals, chrom_size,
                            fragment_length, read_rate)
    peaks = call_peaks(p_values, p_value_cutoff, fragment_length)
    return Interval(intervals.chromosome[:len(peaks)], peaks.start, peaks.end)


def main(filename, genome_file, fragment_length=150, p_value_cutoff=0.001):
    genome = bnp.open(genome_file).read()
    genome_size = genome.size.sum()
    genome_size = 2700000000
    chrom_sizes = GroupedDict((str(name), size) for name, size in zip(genome.name, genome.size))
    intervals = bnp.open(filename, buffer_type=BamIntervalBuffer).read_chunks()
    grouped_intervals = groupby(intervals, "chromosome")
    n_reads = 184961616# bnp.count_entries(filename)
    read_rate = n_reads/genome_size
    return macs2(grouped_intervals, chrom_sizes, fragment_length, read_rate, p_value_cutoff, fragment_length)


with bnp.open("/home/knut/Data/peaks.bed", "w") as outfile:
    for chrom, data in main("/home/knut/Data/ENCFF296OGN.bam", "/home/knut/Data/hg38.chrom.sizes"):
        if len(data):
            outfile.write(data)


#intervals = bnp.open("/home/knut/Data/ENCFF296OGN.bed", buffer_type=Bed6Buffer).read_chunks()


# intervals = bnp.open("/home/knut/Data/ENCFF296OGN.bam", buffer_type=BamIntervalBuffer).read_chunks(chunk_size=524288*32)
# # intervals = alignment_to_interval(reads)
# # first_chunk = next(chunk_lines(iter(intervals), 10000000))
# grouped = groupby(intervals, "chromosome")
# chrom, first_chunk = next(iter(grouped))
# chrom_size = 248956422
# get_p_values(first_chunk, 300, chrom_size)
