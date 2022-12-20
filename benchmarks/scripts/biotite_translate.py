import biotite
import biotite.sequence.io.fasta as fasta
fasta_file = fasta.FastaFile.read(snakemake.input[0])
out_fasta_file = fasta.FastaFile()

for header, dna in fasta_file.items():
    dna = biotite.sequence.NucleotideSequence(dna)
    protein = dna.translate(complete=True)
    fasta.set_sequence(out_fasta_file, protein, header=header)

out_fasta_file.write(snakemake.output[0])
