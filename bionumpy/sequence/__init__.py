from .kmers import get_kmers
from .minimizers import get_minimizers
from .translate import translate_dna_to_protein
from .string_matcher import construct_wildcard_matcher, \
    construct_flexible_len_regex_matchers, \
    construct_fixed_len_regex_matchers
from .position_weight_matrix import get_motif_scores, pwm_from_counts
