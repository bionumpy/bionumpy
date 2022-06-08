"""Top-level package for bionumpy."""

__author__ = """Knut Rand"""
__email__ = "knutdrand@gmail.com"
__version__ = "0.1.0"

from .files import bnp_open as open
from .sequences import Sequence, Sequences, as_sequence_array
from .kmers import KmerEncoding
from .minimizers import Minimizers
from .position_weight_matrix import PositionWeightMatrix
from . import parser
import npstructures as nps
__all__ = ["Sequence", "Sequences", "as_sequence_array",
           "KmerEncoding", "Minimizers", "PositionWeightMatrix"]


def set_backend(cp):
    nps.set_backend(cp)
    parser.wrapper = cp.asarray
    kmers.np = cp
    encodings.np = cp
    delimited_buffers.np = cp
    parser.np = cp
