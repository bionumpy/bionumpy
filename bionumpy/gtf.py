import itertools
import numpy as np
import bionumpy as bnp
from .datatypes import GFFEntry, GFFExonEntry, GFFTranscriptEntry, GFFTranscriptEntry
from .bnpdataclass import bnpdataclass
from .encoded_array import change_encoding
from .encodings import BaseEncoding
from .io.strops import str_equal, split
from . import streamable, EncodedRaggedArray, as_encoded_array
from .sequence import get_reverse_complement
from .datatypes import SequenceEntry

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


def get_stranded_intervals(encoded_array, stranded_intervals):
    pass

@streamable()
def get_transcript_sequences(gtf_entries, reference_sequence):
    if len(gtf_entries) == 0:
        return SequenceEntry.empty()
    reference_sequence = as_encoded_array(reference_sequence, bnp.encodings.alphabet_encoding.ACGTnEncoding)
    exon_entries = get_exons(gtf_entries)
    exon_sequences = reference_sequence[exon_entries.start:exon_entries.stop]
    flat_exon_seqeunece = exon_sequences.ravel()
    groups = itertools.groupby(exon_entries, key=lambda entry: str(entry.transcript_id))
    infos = []
    for transcript_id, entries in groups:
        entries = list(entries)
        strand = str(entries[0].strand)
        seq_length = sum(entry.stop - entry.start for entry in entries)
        infos.append((transcript_id, strand, seq_length))
    names, strands, lengths = zip(*infos)
    transcript_sequences = EncodedRaggedArray(flat_exon_seqeunece, list(lengths))
    transcript_sequences = np.where((as_encoded_array(''.join(strands)) == "-")[:, np.newaxis], get_reverse_complement(transcript_sequences), transcript_sequences)

    # convert transcript encoding to BaseEncoding so that we can create a
    # SequenceEntry.
    return SequenceEntry(list(names), change_encoding(transcript_sequences, BaseEncoding))
