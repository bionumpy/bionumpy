from .kmers import get_kmers, count_kmers
from .minimizers import get_minimizers
from .position_weight_matrix import get_motif_scores, PWM
from .string_matcher import match_string
from .translate import translate_dna_to_protein
from .dna import get_reverse_complement, get_strand_specific_sequences, get_sequences
from .count_encoded import count_encoded, EncodedCounts


__all__ = ["get_kmers", 
           "get_minimizers", 
           "translate_dna_to_protein",
           "get_motif_scores", 
           "get_reverse_complement",
           "get_strand_specific_sequences", 
           "count_encoded",
           "match_string",
           "EncodedCounts",
           ]


def set_backend(lib):
    from . import kmers
    kmers.np = lib
