import bionumpy as bnp
import numpy as np
import matplotlib.pyplot as plt

f = bnp.open("example_data/big.fq.gz")
entries = f.read()
print(entries)
print(len(entries))
plt.hist(entries.sequence.shape.lengths);plt.show()

mean_qualities=entries.quality.mean(axis=-1)

plt.hist(mean_qualities)

plt.show()

gc_mask = (entries.sequence == bnp.as_sequence_array("G")) | (entries.sequence == bnp.as_sequence_array("C"))
gc_content = np.mean(gc_mask, axis=-1)
plt.hist(gc_content); plt.show()
plt.plot(gc_content, mean_qualities, "."); plt.show()

plt.plot(entries.sequence.shape.lengths, mean_qualities, "."); plt.show()

