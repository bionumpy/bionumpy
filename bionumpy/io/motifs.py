from .jaspar import read_jaspar_matrix, read_csv_motif
import numpy as np
from dataclasses import dataclass
from pathlib import PurePath


parsers = {".jaspar": read_jaspar_matrix, ".csv": read_csv_motif}


@dataclass
class Motif:
    alphabet: str
    matrix: np.ndarray


def read_motif(filename):
    path = PurePath(filename)
    suffix = path.suffixes[-1]
    parser = parsers[suffix]
    return parser(filename)
# @ return Motif(alphabet, matrix)
