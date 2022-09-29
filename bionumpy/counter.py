import numpy as np
import dataclasses


@dataclasses.dataclass
class EncodedCounts:
    alphabet: list
    counts: np.ndarray
    
    def __str__(self):
        return "\n".join(f"{c}: {n}" for c, n in zip(self.alphabet, self.counts))


def count_encoded(values):
    alphabet = values.encoding.get_alphabet()
    counts = np.bincount(values, minlength=len(alphabet))
    return EncodedCounts(alphabet, counts)
