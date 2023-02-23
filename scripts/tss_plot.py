import bionumpy as bnp
from bionumpy.bnpdataclass import replace
from bionumpy.io.delimited_buffers import NarrowPeakBuffer
from bionumpy.genomic_data.genomic_intervals import LocationEntry
import matplotlib.pyplot as plt
import numpy as np
from bionumpy.computation_graph import compute
import logging
logging.basicConfig(level=logging.DEBUG)
import typer


def tss_plot(wig_filename: str, chrom_sizes_filename: str, annotation_filename: str):
    genome = bnp.Genome.from_file(chrom_sizes_filename, sort_names=True)
    logging.info('Reading track')
    track = genome.read_track(wig_filename, stream=True)
    logging.info('Reading annotation')
    annotation = genome.read_annotation(annotation_filename)
    transcripts = annotation.transcripts
    tss = transcripts.get_location('start')
    # tss = tss.sorted()
    logging.info('subsetting')
    windows = tss.get_windows(flank=500)
    signals = track[windows]
    print(type(signals))
    logging.info('getting mean')
    mean_signal = signals.mean(axis=0)
    print(type(mean_signal))
    logging.info('computing')
    signal, = compute(mean_signal)
    print(signal)
    signal = signal.to_array()
    plt.plot(np.arange(-500, 500), signal)
    plt.show()


def peak_plot(wig_filename: str, chrom_sizes_filename: str, peak_filename: str):
    genome = bnp.Genome.from_file(chrom_sizes_filename)
    peaks = genome.read_intervals(peak_filename)
    peak_sizes = peaks.stop-peaks.start

    max_len = np.max(peak_sizes)
    mid_points = peaks.get_location('center')
    windows = mid_points.get_windows(flank=(max_len+1)//2)
    #padded_peaks = replace(peaks, start=mid_points-max_len//2, stop=mid_points+max_len//2)
    track = genome.read_track(wig_filename, stream=False)
    signals = track[windows]
    sample_rate = max(max_len//1000, 1)
    signals = signals[:, ::sample_rate]
    signals, = compute(signals)
    args = np.argsort(peak_sizes)
    plt.hist(peak_sizes); plt.show()
    plt.imshow(signals.to_array()[args])
    plt.show()


def summit_plot(bam_filename: str, chrom_sizes_filename: str, peak_filename: str):
    genome = bnp.Genome.from_file(chrom_sizes_filename)
    peaks = bnp.open(peak_filename, buffer_type=NarrowPeakBuffer).read()
    location_entries = LocationEntry(peaks.chromosome, peaks.start+peaks.summit)
    summits = genome.get_locations(location_entries)
    windows = summits.get_windows(flank=200)
    reads = genome.read_intervals(bam_filename, stream=False, stranded=True)
    means = [reads[reads.strand==strand].get_pileup()[windows].mean(axis=0)
             for strand in '+-']

    #pos_reads = reads[reads.strand=='+']
    #neg_reads = reads[reads.strand=='-']
    #pos_pileup = pos_reads.get_pileup()
    #neg_pileup = neg_reads.get_pileup()
    #pos_signals = pos_pileup[windows]
    #neg_signals = neg_pileup[windows]
    pos_mean, neg_mean = compute(*means)
    # pos_mean, neg_mean = compute(pos_signals.mean(axis=0), neg_signals.mean(axis=0))
    plt.plot(np.arange(-200, 200), pos_mean.to_array())
    plt.plot(np.arange(-200, 200), neg_mean.to_array())
    plt.show()


def vcf_plot(wig_filename: str, chrom_sizes_filename: str, vcf_filename: str):
    genome = bnp.Genome.from_file(chrom_sizes_filename)
    variants = genome.get_locations(bnp.open(vcf_filename).read())
    windows = variants.get_windows(flank=100)
    track = genome.read_track(wig_filename, stream=True)
    signals = track[windows]
    mean_signal = signals.mean(axis=0)
    signal, = compute(mean_signal)
    signal = signal.to_array()
    plt.plot(np.arange(-500, 500), signal)
    plt.show()


def main(wig_filename: str, chrom_sizes_filename: str, filename: str):
    if wig_filename.endswith('.bam')or wig_filename.endswith('reads.bed'):
        func = summit_plot
    elif filename.endswith('vcf.gz'):
        func = vcf_plot
    elif filename.endswith('bed.gz'):
        func = peak_plot
    else:
        func = tss_plot
    func(wig_filename, chrom_sizes_filename, filename)


# tss_plot(*('/home/knut/Data/out.wig /home/knut/Data/hg38.chrom.sizes /home/knut/Data/gencode.v43.annotation.gff3.gz'.split()))


#main(*('/home/knut/Data/out.wig /home/knut/Data/hg38.chrom.sizes /home/knut/Data/ENCFF266FSE.bed.gz'.split()))
# main(*'example_data/CTCF_chr21-22.wig example_data/chr21-22.chrom.sizes example_data/chr21a22.gtf'.split())
# main(*'example_data/CTCFpvalues_chr21-22.wig example_data/chr21-22.chrom.sizes example_data/ctcf_chr21-22.bed.gz'.split())
main(*'example_data/ctcf_chr21-22_reads.bed example_data/chr21-22.chrom.sizes example_data/ctcf_chr21-22.bed.gz'.split())
# if __name__ == '__main__':
#     typer.run(main)
