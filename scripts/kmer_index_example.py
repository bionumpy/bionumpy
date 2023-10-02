from bionumpy.io.delimited_buffers import DelimitedBuffer, get_bufferclass_for_datatype
from bionumpy.bnpdataclass import bnpdataclass
import bionumpy as bnp
from bionumpy.sequence.indexing import KmerLookup


@bnpdataclass
class Olga:
    dna: bnp.DNAEncoding
    amino_acid: bnp.AminoAcidEncoding


class OlgaBuffer(DelimitedBuffer):
    dataclass = Olga


def test():
    olga_sequence_data = bnp.open(filename="example_data/airr.tsv", buffer_type=OlgaBuffer).read()
    dna_3mer_lookup = KmerLookup.create_lookup(olga_sequence_data.dna, k=3)
    tgc_sequences = dna_3mer_lookup.get_sequences("TGC")
    print(tgc_sequences)
