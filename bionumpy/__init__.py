"""Top-level package for bionumpy."""

__author__ = """Knut Rand"""
__email__ = "knutdrand@gmail.com"
__version__ = "0.1.0"

from .files import bnp_open as open
from .files import count_entries
from .sequences import Sequence, as_sequence_array, as_encoded_sequence_array
from .kmers import KmerEncoding
from .minimizers import Minimizers
from .position_weight_matrix import PositionWeightMatrix
from .counter import count_encoded
from . import parser
from .npdataclassstream import mean, bincount, histogram, streamable
import npstructures as nps
from .file_buffers import (TwoLineFastaBuffer, FastQBuffer)
from .delimited_buffers import (BedBuffer, VCFBuffer, VCFMatrixBuffer,
                                GfaSequenceBuffer, get_bufferclass_for_datatype)
from .datatypes import SAMEntry, GFFEntry, Bed6
from .multiline_buffer import MultiLineFastaBuffer
from .encodings.alphabet_encoding import (DNAEncoding, RNAENcoding, AminoAcidEncoding,
                                          DNAArray, RNAArray, AminoAcidArray)

SAMBuffer = get_bufferclass_for_datatype(SAMEntry)
GFFBuffer = get_bufferclass_for_datatype(GFFEntry)
Bed6Buffer = get_bufferclass_for_datatype(Bed6)

__all__ = ["Sequence", "as_sequence_array", "as_encoded_sequence_array",
           "KmerEncoding", "Minimizers", "PositionWeightMatrix", "mean",
           "bincount", "streamable", "histogram", "count_entries",
           "BedBuffer", "VCFBuffer", "VCFMatrixBuffer", "GfaSequenceBuffer",
           "TwoLineFastaBuffer", "FastQBuffer",
           "SAMBuffer", "GFFBuffer", "Bed6Buffer", "MultiLineFastaBuffer",
           "count_encoded", "DNAEncoding", "RNAENcoding", "AminoAcidEncoding"]


def set_backend(cp):
    nps.set_backend(cp)
    parser.wrapper = cp.asarray
    kmers.np = cp
    encodings.np = cp
    delimited_buffers.np = cp
    parser.np = cp
