from .jaspar import read_jaspar_matrix
import numpy as np
from dataclasses import dataclass
from pathlib import PurePath


parsers = {".jaspar": read_jaspar_matrix}


@dataclass
class Motif:
    alphabet: str
    matrix: np.ndarray


def read_motif(filename):
    path = PurePath(filename)
    suffix = path.suffixes[-1]
    parser = parsers[suffix]
    alphabet, matrix = parser(filename)
    return Motif(alphabet, matrix)
