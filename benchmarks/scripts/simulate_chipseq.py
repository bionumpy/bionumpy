from pyjaspar import jaspardb
import bionumpy as bnp
from bionumpy.sequence.position_weight_matrix import PWM
from bionumpy.datatypes import ChromosomeSize, SequenceEntry
from bionumpy.simulate.chipseq import simulate_chip_seq_reads, simulate_sequence, ChipSeqSimulationSettings
jaspar_object = jaspardb(release="JASPAR2020")
ctcf_motif = jaspar_object.fetch_motifs_by_name('CTCF')[0]
motif = PWM.from_dict(ctcf_motif.pwm)
chromosome_sizes = {"chr1": 10000,
                    "chr2": 20000,
                    "chr3": 30000}
chrom_sizes = ChromosomeSize(list(chromosome_sizes.keys()),
                             list(chromosome_sizes.values()))

# motif = read_motif("example_data/MA0080.1.jaspar")
settings = ChipSeqSimulationSettings(motif)
sequences = {name: simulate_sequence("acgt", size) for name, size in chromosome_sizes.items()}
multistream = bnp.MultiStream(chromosome_sizes, sequences=sequences)
reads = simulate_chip_seq_reads(multistream.sequences, settings, multistream.sequence_names)
with bnp.open(snakemake.output[2], 'w') as f:
    f.write(SequenceEntry(list(sequences.keys()), list(sequences.values())))
with bnp.open(snakemake.output[1], 'w', buffer_type=bnp.Bed6Buffer) as f:
    f.write(reads)

with bnp.open(snakemake.output[0], "w") as f:
    f.write(chrom_sizes)
        
