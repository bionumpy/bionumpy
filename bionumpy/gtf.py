import itertools
import numpy as np
import bionumpy as bnp
from .datatypes import GFFEntry
from .bnpdataclass import bnpdataclass
from .io.strops import str_equal, split
from . import streamable, EncodedRaggedArray, as_encoded_array
from .dna import reverse_compliment
from .datatypes import SequenceEntry


@bnpdataclass
class GFFGeneEntry(GFFEntry):
    gene_id: str


@bnpdataclass
class GFFTranscriptEntry(GFFGeneEntry):
    transcript_id: str


@bnpdataclass
class GFFExonEntry(GFFTranscriptEntry):
    exon_id: str

def get_attributes(gtf_entries, attribute_names):
    gtf_entries.atributes[:, -1] = " "
    all_features = split(gtf_entries.atributes.ravel(), " ")
    keys = all_features[:-1:2]
    values = all_features[1::2, 1:-2]
    return {name: values[str_equal(keys, name)] for name in attribute_names}


@streamable()
def get_genes(gtf_entries):
    genes = gtf_entries[str_equal(gtf_entries.feature_type, "gene")]
    attributes = get_attributes(genes, ["gene_id"])
    for key, val in attributes.items():
        assert len(val) == len(genes), (len(val), len(genes))
    return GFFGeneEntry(*genes.shallow_tuple(), **attributes)


@streamable()
def get_transcripts(gtf_entries):
    genes = gtf_entries[str_equal(gtf_entries.feature_type, "transcript")]
    attributes = get_attributes(genes, ["transcript_id", "gene_id"])
    return GFFTranscriptEntry(*genes.shallow_tuple(), **attributes)


@streamable()
def get_exons(gtf_entries):
    genes = gtf_entries[str_equal(gtf_entries.feature_type, "exon")]
    attributes = get_attributes(genes, ["transcript_id", "gene_id", "exon_id"])
    return GFFExonEntry(*genes.shallow_tuple(), **attributes)


@streamable()
def get_transcript_sequences(gtf_entries, reference_sequence):
    reference_sequence = as_encoded_array(reference_sequence, bnp.encodings.alphabet_encoding.ACGTnEncoding)
    exon_entries = get_exons(gtf_entries)
    exon_sequences = reference_sequence[exon_entries.start:exon_entries.stop]
    groups = itertools.groupby(zip(exon_entries, exon_sequences), key=lambda pair: str(pair[0].transcript_id))
    sequences = {}
    for transcript_id, pairs in groups:
        pairs = list(pairs)
        strand = pairs[0][0].strand
        sequence = np.concatenate([pair[1] for pair in pairs])
        if strand == "-":
            sequence = reverse_compliment(sequence)
        sequence = EncodedRaggedArray(sequence, [len(sequence)])
        sequences[transcript_id] = sequence
    return SequenceEntry(list(sequences.keys()), np.concatenate(list(sequences.values())))
