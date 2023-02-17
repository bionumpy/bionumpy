import typer
import bionumpy as bnp
import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def main(fasta_filename: str, annotation_filename: str):
    genome = bnp.Genome.from_file(fasta_filename)
    logging.info('Reading annotation')
    annotation = genome.read_annotation(annotation_filename)
    logging.info('Reading sequence')
    reference_sequence = genome.read_sequence()
    logging.info('Calculating masks')
    transcription_mask = annotation.transcripts.get_mask()
    exon_mask = annotation.exons.get_mask()
    intron_mask = transcription_mask & ~exon_mask
    logging.info('Extracting sequences')
    exon_sequence = reference_sequence[exon_mask]
    intron_sequence = reference_sequence[intron_mask]
    print(bnp.count_encoded(exon_sequence))
    print(bnp.count_encoded(intron_sequence))


if __name__ == '__main__':
    typer.run(main)
