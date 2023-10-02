from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

with open(snakemake.output[0],'w') as aa_fa:
    for dna_record in SeqIO.parse(snakemake.input[0],'fasta'):
        new_record = SeqRecord(
            dna_record.seq.reverse_complement(),
            id=dna_record.id, description="")
        SeqIO.write(new_record, aa_fa,'fasta-2line')
