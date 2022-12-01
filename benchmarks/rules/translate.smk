import biotite.sequence

import bionumpy as bnp


rule translate_bionumpy:
    input:
        "results/dna_sequences/{filename}.fa"
    output:
        "results/bionumpy/translate/{filename}.fa"
    benchmark:
        "benchmarks/translate/bionumpy/{filename}.txt"
    run:
        from bionumpy.sequence import translate_dna_to_protein
        input_stream = bnp.open(input[0]).read_chunks()
        output_stream = bnp.open(output[0], "w")
        output_stream.write(translate_dna_to_protein(input_stream))
        output_stream.close()

rule translate_biopython:
    input:
        "results/dna_sequences/{name}.fa"
    output:
        "results/biopython/translate/{name}.fa"
    benchmark:
        "benchmarks/translate/biopython/{name}.txt"
    run:
        from Bio import SeqIO
        from Bio.SeqRecord import SeqRecord
        with open(output[0], 'w') as aa_fa:
            for dna_record in SeqIO.parse(input[0], 'fasta'):
                new_record = SeqRecord(
                    dna_record.seq.translate(),
                    id=dna_record.id)
                SeqIO.write(new_record, aa_fa, 'fasta')


rule translate_biostrings:
    input:
        "results/dna_sequences/{name}.fa"
    output:
        "results/biostrings/translate/{name}.fa"
    benchmark:
        "benchmarks/translate/biostrings/{name}.txt"
    script:
        "scripts/reverse_complement_biostrings.R"


def translate_dna_to_protein(seq):
    # source: https://www.geeksforgeeks.org/dna-protein-python-3/
    table = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
        'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
    }
    protein = ""
    if len(seq) % 3 == 0:
        for i in range(0,len(seq),3):
            codon = seq[i:i + 3]
            protein += table[codon]
    return protein


rule translate_python:
    input:
        "results/dna_sequences/{name}.fa"
    output:
        "results/python/translate/{name}.fa"
    benchmark:
        "benchmarks/translate/python/{name}.txt"
    run:
        with open(input[0]) as infile:
            with open(output[0], "w") as outfile:
                for line in infile:
                    if line.startswith(">"):
                        outfile.write(line)
                    else:
                        outfile.write(
                            translate_dna_to_protein(line.strip()) + "\n")



rule translate_biotite:
    input:
        "results/dna_sequences/{name}.fa"
    output:
        "results/biotite/translate/{name}.fa"
    benchmark:
        "benchmarks/translate/biotite/{name}.txt"

    run:
        import biotite.sequence.io.fasta as fasta
        fasta_file = fasta.FastaFile.read(input[0])
        out_fasta_file = fasta.FastaFile()

        for header, dna in fasta_file.items():
            dna = biotite.sequence.NucleotideSequence(dna)
            protein = dna.translate(complete=True)
            fasta.set_sequence(out_fasta_file, protein, header=header)

        out_fasta_file.write(output[0])
