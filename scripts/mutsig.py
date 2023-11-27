import bionumpy as bnp
from bionumpy.variants import count_mutation_types_genomic
from bionumpy.io.matrix_dump import matrix_to_csv
from bionumpy import Genome

import logging
logging.basicConfig(level="INFO")


def main(vcf_filename: str, fasta_filename: str, out_filename: str = None, has_numeric_chromosomes=True):
    genome = Genome.from_file(fasta_filename)
    sequence = genome.read_sequence()
    counts = None
    for chunk in bnp.open(vcf_filename).read_chunks():
        variants = genome.get_locations(chunk, has_numeric_chromosomes=has_numeric_chromosomes)
        if counts is None:
            counts = count_mutation_types_genomic(variants, sequence, genotyped=False)
        else:
            counts += count_mutation_types_genomic(variants, sequence, genotyped=False)

    if out_filename is not None:
        output = matrix_to_csv(counts.counts, header=counts.alphabet)
        open(out_filename, "ab").write(bytes(output.raw()))
    return counts


if __name__ == "__main__":
    import sys
    import os
    in_folder = sys.argv[1]
    fasta_file = sys.argv[2]
    out_file = sys.argv[3]
    f = open(out_file, 'w')
    written_header = False
    for filename in os.listdir(in_folder):
        if not (filename.endswith(".vcf") or filename.endswith(".vcf.gz")):
            continue
        logging.info('Processing %s', filename)
        counts = main(os.path.join(in_folder, filename), fasta_file)
        if not written_header:
            f.write('filename' + ',' + ','.join(counts.alphabet) + '\n')
            written_header = True
        f.write(filename + ',' + ','.join(map(str, counts.counts)) + '\n')


