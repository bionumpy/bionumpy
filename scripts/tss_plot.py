import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import plotly.express as px
import bionumpy as bnp
import typer
plot=True


def tss_plot(wig_filename: str, chrom_sizes_filename: str, annotation_filename: str, plot=True):
    # Read genome and transcripts
    genome = bnp.Genome.from_file(chrom_sizes_filename, sort_names=True) # The wig file is alphbetically sorted
    annotation = genome.read_annotation(annotation_filename)
    transcripts = annotation.transcripts

    # Get transcript start locations and make windows around them
    tss = transcripts.get_location('start').sorted() # Make sure the transcripts are sorted alphabetically
    windows = tss.get_windows(flank=500)

    # Get mean read pileup within these windows and plot
    track = genome.read_track(wig_filename, stream=True)
    signals = track[windows]
    mean_signal = signals.mean(axis=0)
    signal = bnp.compute(mean_signal)  # Compute the actual value
    if plot:
        return px.line(x=np.arange(-500, 500), y=signal.to_array(),
                       title="Read pileup relative to TSS start",
                       labels={"x": "Position relative to TSS start", "y": "Mean read pileup"})


def cpg_plot(fasta_filename: str, annotation_filename: str, plot: bool=True):
    genome = bnp.Genome.from_file(fasta_filename, sort_names=False)
    reference_sequence = genome.read_sequence()
    annotation = genome.read_annotation(annotation_filename)
    transcripts = annotation.transcripts
    # Get transcript start locations and make windows around them
    tss = transcripts.get_location('start')
    tss = tss.sorted() # Make sure the transcripts are sorted alphabetically
    flank = 200
    windows = tss.get_windows(flank=flank)
    cpg_proportion = bnp.match_string(reference_sequence[windows], 'CG').mean(axis=0)
    if plot:
        return px.line(x=np.arange(-flank, flank-1), y=cpg_proportion,
                       title="Read pileup relative to TSS start",
                       labels={"x": "Position relative to TSS start", "y": "Mean read pileup"})


def cpg_tss_plot(wig_filename: str, fasta_filename: str, annotation_filename: str, plot=True):
    # Read genome and transcripts
    genome = bnp.Genome.from_file(fasta_filename, sort_names=False) # The wig file is alphbetically sorted
    reference_sequence = genome.read_sequence()
    annotation = genome.read_annotation(annotation_filename)
    transcripts = annotation.transcripts
    # Get transcript start locations and make windows around them
    tss = transcripts.get_location('start')
    tss = tss.sorted() # Make sure the transcripts are sorted alphabetically
    windows = tss.get_windows(flank=500)
    assert windows.is_stranded()
    cpg_proportion = bnp.match_string(reference_sequence[windows], 'CG').mean(axis=-1)
    has_cpg = np.log2(cpg_proportion) > -5
    # px.histogram(np.log2(cpg_proportion)).show()
    # Get mean read pileup within these windows and plot
    track = genome.read_track(wig_filename, stream=True)
    signals = bnp.compute(track[windows])
    mean_signals = {'cpg': signals[has_cpg].mean(axis=0),
                    'non-cpg': signals[~has_cpg].mean(axis=0)}
    # mean_signals = bnp.compute(mean_signals)  # Compute the actual value
    if plot:
        return go.Figure(
            [go.Scatter(x=np.arange(-500, 500), y=signal.to_array(), name=name)
             for name, signal in mean_signals.items()],
            layout={'title': "Read pileup relative to TSS start",
                    'xaxis_title': "Position relative to TSS start",
                    'yaxis_title': 'Read coverage'})
        
        # px.line(x=np.arange(-500, 500), y=signal.to_array(),
        #         title=,
        #         labels={"x": , "y": "Mean read pileup"}).show()


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
    signals = bnp.compute(signals)
    args = np.argsort(peak_sizes)
    plt.imshow(signals.to_array()[args])
    plt.show()


def summit_plot(bam_filename: str, chrom_sizes_filename: str, peak_filename: str, plot=True):

    # Read genome and peaks
    genome = bnp.Genome.from_file(chrom_sizes_filename).with_ignored_added(['chrEBV'])
    peaks = bnp.open(peak_filename, buffer_type=bnp.NarrowPeakBuffer).read()
    location_entries = bnp.LocationEntry(peaks.chromosome, peaks.start+peaks.summit)
    # Create locations of peaks summits
    summits = genome.get_locations(location_entries).sorted()

    # Create windows around summits and extract read pileup
    windows = summits.get_windows(flank=200)
    reads = genome.read_intervals(bam_filename, stream=False, stranded=True)

    # Get mean pileup for reads with negative and positive strand
    signals_dict = {strand: reads[reads.strand == strand].get_pileup()[windows].mean(axis=0)
                    for strand in '+-'}
    signals_dict = bnp.compute(signals_dict)
    if plot:
        return go.Figure(
            [go.Scatter(x=np.arange(-200, 200), y=signal.to_array(), name=f'{strand} Strand')
             for strand, signal in signals_dict.items()],
            layout={'title': 'Summit plot',
                    'xaxis_title': 'Distance from peak summit',
                    'yaxis_title': 'Read coverage'})
        # fig.add_trace(go.Scatter(x=np.arange(-200, 200), y=pos_mean.to_array(), name='Positive Strand'))
        #fig.add_trace(go.Scatter(x=np.arange(-200, 200), y=neg_mean.to_array(), name='Negative Strand'))
        # fig.show()


def vcf_plot(wig_filename: str, chrom_sizes_filename: str, vcf_filename: str, plot=True):
    # Read genome and variants
    genome = bnp.Genome.from_file(chrom_sizes_filename)
    variants = genome.read_locations(vcf_filename, has_numeric_chromosomes=True)

    # Get windows around variants and get read pileup in these windows
    flank = 100
    windows = variants.get_windows(flank=flank)
    reads = genome.read_intervals(wig_filename, stream=True, stranded=True)
    track = reads.get_pileup()
    signals = track[windows]

    # Get mean signal inside these windows and plot
    mean_signal = signals.mean(axis=0)
    signal = bnp.compute(mean_signal)
    signal = signal.to_array()

    if plot:
        return px.line(x=np.arange(-flank, flank), y=signal,
                       title="Read pileup relative to common variants",
                       labels={"x": "Position relative to variant location", "y": "Mean read pileup"})


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


def test(plot=False):
    tss_plot("example_data/CTCF_chr21-22.wig.gz", "example_data/chr21-22.chrom.sizes", "example_data/chr21a22.gtf", plot=plot)
    summit_plot("example_data/ctcf_chr21-22.bam", "example_data/chr21-22.chrom.sizes", "example_data/ctcf_chr21-22.bed.gz", plot=plot)
    vcf_plot('example_data/ctcf_chr21-22.bam', 'example_data/chr21-22.chrom.sizes', 'example_data/1000Genomes_chr21-22.vcf.gz', plot=plot)
# tss_plot(*('/home/knut/Data/out.wig /home/knut/Data/hg38.chrom.sizes /home/knut/Data/gencode.v43.annotation.gff3.gz'.split()))


def write_images():
    #fig = tss_plot(*('/home/knut/Downloads/ENCFF994MIH.wig /home/knut/Data/hg38.chrom.sizes /home/knut/Data/gencode.v43.annotation.gff3.gz'.split()), True)
    # fig.write_image('docs_source/figures/tss_plot_ENCFF994MIH.png')
    fig = summit_plot("example_data/ctcf_chr21-22.bam", "example_data/chr21-22.chrom.sizes", "example_data/ctcf_chr21-22.bed.gz", plot=plot)
    fig.write_image('docs_source/figures/summit_plot.png')
    fig = vcf_plot('example_data/ctcf_chr21-22.bam', 'example_data/chr21-22.chrom.sizes', 'example_data/1000Genomes_chr21-22.vcf.gz', plot=plot)
    fig.write_image('docs_source/figures/vcf_plot.png')


write_images()
# main(*('/home/knut/Data/out.wig /home/knut/Data/hg38.chrom.sizes /home/knut/Data/ENCFF266FSE.bed.gz'.split()))
#main(*'example_data/CTCF_chr21-22.wig example_data/chr21-22.chrom.sizes example_data/chr21a22.gtf'.split())
# main(*'example_data/CTCFpvalues_chr21-22.wig example_data/chr21-22.chrom.sizes example_data/ctcf_chr21-22.bed.gz'.split())
#main(*'example_data/ctcf_chr21-22_reads.bed example_data/chr21-22.chrom.sizes example_data/ctcf_chr21-22.bed.gz'.split())
# vcf_plot('example_data/ctcf_chr21-22.bam', 'example_data/chr21-22.chrom.sizes', 'example_data/1000Genomes_chr21-22.vcf.gz')
# cpg_plot('/home/knut/Sources/bionumpy/example_data/sacCer3.fa',
#          '/home/knut/Sources/bionumpy/example_data/sacCer3.ensGene.gtf.gz', plot=False)
#cpg_plot('/home/knut/Data/hg38.fa', '/home/knut/Data/gencode.v41.annotation.gtf.gz',
#          plot=False)

# if __name__ == '__main__':

#     typer.run(cpg_plot)
