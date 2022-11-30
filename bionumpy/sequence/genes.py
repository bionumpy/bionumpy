import itertools
import numpy as np
from ..streams import streamable
from ..sequence import get_reverse_complement
from ..datatypes import SequenceEntry
from ..encoded_array import as_encoded_array, EncodedRaggedArray, change_encoding
from ..encodings import BaseEncoding
from ..encodings.alphabet_encoding import ACGTnEncoding



@streamable()
def get_transcript_sequences(gtf_entries, reference_sequence):
    if len(gtf_entries) == 0:
        return SequenceEntry.empty()
    reference_sequence = as_encoded_array(reference_sequence, ACGTnEncoding)
    exon_entries = gtf_entries.get_exons()
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
    transcript_sequences = np.where((as_encoded_array(''.join(strands)) == "-")[:, np.newaxis], 
                                    get_reverse_complement(transcript_sequences), transcript_sequences)
    return SequenceEntry(list(names), change_encoding(transcript_sequences, BaseEncoding))
