import numpy as np

from bionumpy import DNAEncoding
from bionumpy.simulate import simulate_sequences
from bionumpy.sequence.string_matcher import StringMatcher

named_seqs = simulate_sequences('ACGT', {f's{i}':10+i for i in range(20)})
print(named_seqs)
seqs = named_seqs.sequence
gc_content_per_seq = np.mean((seqs=='C')|(seqs=='G'), axis=1)
print(gc_content_per_seq)
named_seqs = named_seqs.add_fields({'gc':gc_content_per_seq}, {'gc':float})
print(named_seqs)
ac_hits = StringMatcher("AC",DNAEncoding).rolling_window(seqs)
ac_hit_sums = np.sum(ac_hits,axis=1)
print(ac_hit_sums)
named_seqs = named_seqs.add_fields({'ac_hits':ac_hit_sums}, {'ac_hits':int})
print(named_seqs)