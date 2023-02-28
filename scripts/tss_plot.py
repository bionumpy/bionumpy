import numpy as np
import matplotlib.pyplot as plt
import typer
import bionumpy as bnp


def tss_plot(wig_filename: str, chrom_sizes_filename: str, annotation_filename: str):
    genome = bnp.Genome.from_file(chrom_sizes_filename, sort_names=True)
    track = genome.read_track(wig_filename, stream=True)
    annotation = genome.read_annotation(annotation_filename)
    transcripts = annotation.transcripts
    tss = transcripts.get_location('start').sorted()
    windows = tss.get_windows(flank=500)
    signals = track[windows]
    mean_signal = signals.mean(axis=0)
    signal, = bnp.compute(mean_signal)
    signal = signal.to_array()
    plt.plot(np.arange(-500, 500), signal)
    plt.show()


def peak_plot(wig_filename: str, chrom_sizes_filename: str, peak_filename: str):
    genome = bnp.Genome.from_file(chrom_sizes_filename, sort_names=True)
    peaks = genome.read_intervals(peak_filename).sorted()
    peak_sizes = peaks.stop-peaks.start
    max_len = np.max(peak_sizes)
    mid_points = peaks.get_location('center')
    windows = mid_points.get_windows(flank=(max_len+1)//2)
    track = genome.read_track(wig_filename, stream=True)
    signals = track[windows]
    sample_rate = max(max_len//1000, 1)
    signals = signals[:, ::sample_rate]
    signals, = bnp.compute(signals)
    args = np.argsort(peak_sizes)
    plt.imshow(signals.to_array()[args])
    plt.show()


def summit_plot(bam_filename: str, chrom_sizes_filename: str, peak_filename: str):
    genome = bnp.Genome.from_file(chrom_sizes_filename).with_ignored_added(['chrEBV'])
    peaks = bnp.open(peak_filename, buffer_type=bnp.NarrowPeakBuffer).read()
    location_entries = bnp.LocationEntry(peaks.chromosome, peaks.start+peaks.summit)
    summits = genome.get_locations(location_entries).sorted()
    windows = summits.get_windows(flank=200)
    reads = genome.read_intervals(bam_filename, stream=True, stranded=True)
    means = [reads[reads.strand == strand].get_pileup()[windows].mean(axis=0)
             for strand in '+-']
    pos_mean, neg_mean = bnp.compute(*means)
    plt.plot(np.arange(-200, 200), pos_mean.to_array())
    plt.plot(np.arange(-200, 200), neg_mean.to_array())
    plt.show()


def vcf_plot(wig_filename: str, chrom_sizes_filename: str, vcf_filename: str):
    flank = 100
    genome = bnp.Genome.from_file(chrom_sizes_filename)
    variants = genome.read_locations(vcf_filename, has_numeric_chromosomes=True)
    windows = variants.get_windows(flank=flank)
    reads = genome.read_intervals(wig_filename, stream=True, stranded=True)
    track = reads.get_pileup()
    signals = track[windows]
    mean_signal = signals.sum(axis=0)
    signal, = bnp.compute(mean_signal)
    signal = signal.to_array()
    plt.plot(np.arange(-flank, flank), signal)
    plt.show()


def main(wig_filename: str, chrom_sizes_filename: str, filename: str):
    if wig_filename.endswith('.bam') or wig_filename.endswith('reads.bed'):
        func = summit_plot
    elif filename.endswith('vcf.gz'):
        func = vcf_plot
    elif filename.endswith('bed.gz'):
        func = peak_plot
    else:
        func = tss_plot
    func(wig_filename, chrom_sizes_filename, filename)


# tss_plot(*('/home/knut/Data/out.wig /home/knut/Data/hg38.chrom.sizes /home/knut/Data/gencode.v43.annotation.gff3.gz'.split()))


# main(*('/home/knut/Data/out.wig /home/knut/Data/hg38.chrom.sizes /home/knut/Data/ENCFF266FSE.bed.gz'.split()))
main(*'example_data/CTCF_chr21-22.wig example_data/chr21-22.chrom.sizes example_data/chr21a22.gtf'.split())
# main(*'example_data/CTCFpvalues_chr21-22.wig example_data/chr21-22.chrom.sizes example_data/ctcf_chr21-22.bed.gz'.split())
#main(*'example_data/ctcf_chr21-22_reads.bed example_data/chr21-22.chrom.sizes example_data/ctcf_chr21-22.bed.gz'.split())
# vcf_plot('example_data/ctcf_chr21-22.bam', 'example_data/chr21-22.chrom.sizes', 'example_data/1000Genomes_chr21-22.vcf.gz')
if __name__ == '__main__':
    typer.run(main)
