import typer
import bionumpy as bnp
import plotly.graph_objects as go


def bar(counter):
    return go.Figure(go.Bar(
        x=counter.counts,
        y=counter.alphabet,
        orientation='h'))


def main(fasta_filename: str, annotation_filename: str):
    genome = bnp.Genome.from_file(fasta_filename)
    annotation = genome.read_annotation(annotation_filename)
    reference_sequence = genome.read_sequence()
    transcription_mask = annotation.transcripts.get_mask()
    exon_mask = annotation.exons.get_mask()
    intron_mask = transcription_mask & ~exon_mask
    counts = bnp.count_encoded(reference_sequence[exon_mask])
    print(counts)
    intron_sequence = reference_sequence[intron_mask]
    print(intron_sequence)
    bar(bnp.count_encoded(intron_sequence)).show()


if __name__ == '__main__':
    typer.run(main)
