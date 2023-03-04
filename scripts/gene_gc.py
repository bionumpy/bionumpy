import typer
import bionumpy as bnp
import plotly.graph_objects as go


def plot_bars(**counter_dict):
    fig = go.Figure()
    for name, counter in counter_dict.items():
        fig.add_trace(go.Bar(x=counter.alphabet,
                             y=counter.proportions,
                             name=name))
    return fig


def main(fasta_filename: str, annotation_filename: str):
    genome = bnp.Genome.from_file(fasta_filename)
    annotation = genome.read_annotation(annotation_filename)
    reference_sequence = genome.read_sequence()
    transcript_mask = annotation.transcripts.get_mask()
    exon_mask = annotation.exons.get_mask()
    intron_mask = transcript_mask & ~exon_mask
    exon_counts = bnp.count_encoded(reference_sequence[exon_mask])
    intron_counts = bnp.count_encoded(reference_sequence[intron_mask])
    if plot:
        plot_bars(exon=exon_counts, intron=intron_counts).show()


if __name__ == '__main__':
    typer.run(main)
