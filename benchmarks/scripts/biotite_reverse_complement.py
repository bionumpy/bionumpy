import biotite
import biotite.sequence.io.fasta as fasta

fasta_file = fasta.FastaFile.read(snakemake.input[0])
out_fasta_file = fasta.FastaFile(chars_per_line=100000000)

for header, dna in fasta_file.items():
    dna = biotite.sequence.NucleotideSequence(dna)
    revcomp = dna.reverse().complement()
    fasta.set_sequence(out_fasta_file, revcomp, header=header)

out_fasta_file.write(snakemake.output[0])
