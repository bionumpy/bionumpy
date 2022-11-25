import bionumpy as bnp
import logging
from bionumpy.gtf import get_transcript_sequences
from bionumpy.cli import run_as_commandline


def extract_transcriptome(reference_file: str, gtf_file: str, out_file: str):
    logging.basicConfig(level=logging.DEBUG)
    reference = bnp.open_indexed(reference_file)
    annotation = bnp.open(gtf_file).read_chunks()
    multistream = bnp.streams.multistream.MultiStream(reference.get_contig_lengths(), entries= annotation, sequence=reference)
    result = get_transcript_sequences(multistream.entries, multistream.sequence)
    with bnp.open(out_file, "w") as outfile:
        outfile.write(result)


def test():
    pass

if __name__ == '__main__':
    run_as_commandline(extract_transcriptome)
    #extract_transcriptome(reference_file="/Users/kanduric/Documents/Projects/christin_transcriptomics/GRCh38.p13.genome.fa",
    #gtf_file="/Users/kanduric/Documents/Projects/christin_transcriptomics/gencode.v38.chr_patch_hapl_scaff.annotation.gtf",
    # out_file="/Users/kanduric/Documents/Projects/christin_transcriptomics/bnp_extracted_transcriptome.fa")
