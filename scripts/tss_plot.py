import bionumpy as bnp
from bionumpy.bnpdataclass import replace
import matplotlib.pyplot as plt
import numpy as np
from bionumpy.computation_graph import compute
import logging
logging.basicConfig(level=logging.DEBUG)
# import typer


def tss_plot(wig_filename: str, chrom_sizes_filename: str, annotation_filename: str):
    genome = bnp.Genome.from_file(chrom_sizes_filename, sort_names=True)
    track = genome.read_track(wig_filename, stream=True)
    annotation = genome.read_annotation(annotation_filename)
    transcripts = annotation.transcripts
    tss = transcripts.get_location('start')
    # tss = tss.sorted()
    windows = tss.get_windows(flank=500)

    signals = track[windows]
    mean_signal = signals.mean(axis=0)
    signal, = compute(mean_signal)
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
    track = genome.read_track(wig_filename, stream=True)
    signals = track[windows]
    sample_rate = max(max_len//1000, 1)
    signals = signals[:, ::sample_rate]
    signals, = compute(signals)
    args = np.argsort(peak_sizes)
    plt.imshow(signals.to_array()[args])
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
    if filename.endswith('vcf.gz'):
        func = vcf_plot
    elif filename.endswith('bed.gz'):
        func = peak_plot
    else:
        func = tss_plot
    func(wig_filename, chrom_sizes_filename, filename)


tss_plot(*('/home/knut/Data/out.wig /home/knut/Data/hg38.chrom.sizes /home/knut/Data/gencode.v43.annotation.gff3.gz'.split()))


# main(*('/home/knut/Data/out.wig /home/knut/Data/hg38.chrom.sizes /home/knut/Data/ENCFF266FSE.bed.gz'.split()))


# if __name__ == '__main__':
#    typer.run(main)
