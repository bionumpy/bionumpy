import bionumpy as bnp
import matplotlib.pyplot as plt
from bionumpy.genomic_data.genomic_intervals import GenomicLocation
import numpy as np
from bionumpy.computation_graph import compute
import typer


def tss_plot(wig_filename: str, chrom_sizes_filename: str, annotation_filename: str):
    genome = bnp.Genome.from_file(chrom_sizes_filename)
    annotation = genome.read_annotation(annotation_filename)
    transcripts = annotation.transcripts
    tss = transcripts.get_location('start')
    windows = tss.get_windows(flank=500)
    track = genome.read_track(wig_filename, stream=True)
    signals = track[windows]
    mean_signal = signals.mean(axis=0)
    signal, = compute(mean_signal)
    signal = signal.to_array()
    # plt.plot(np.arange(-500, 500), signal)
    # plt.show()


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
    else:
        func = tss_plot
    func(wig_filename, chrom_sizes_filename, filename)


tss_plot(*('/home/knut/Data/out.wig /home/knut/Data/hg38.chrom.sizes /home/knut/Data/gencode.v43.annotation.gff3.gz'.split()))
if __name__ == '__main__':
    typer.run(main)
