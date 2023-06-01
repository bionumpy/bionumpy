"""Top-level package for bionumpy."""

__author__ = """Knut Rand"""
__email__ = "knutdrand@gmail.com"
__version__ = '0.2.25'

import npstructures as nps

from .io import (count_entries, open_indexed, MultiLineFastaBuffer, bnp_open,
                 TwoLineFastaBuffer, FastQBuffer, Bed6Buffer, NarrowPeakBuffer,
                 BedBuffer, VCFBuffer, PhasedVCFMatrixBuffer, VCFMatrixBuffer,
                 GfaSequenceBuffer, get_bufferclass_for_datatype)
from .encodings.alphabet_encoding import (DNAEncoding, RNAENcoding, AminoAcidEncoding)
from .encoded_array import EncodedArray, EncodedRaggedArray, as_encoded_array, OneToOneEncoding, BaseEncoding, change_encoding
from .sequence import (get_kmers, get_minimizers, get_motif_scores, count_encoded, match_string)
from .streams import mean, bincount, histogram, streamable, quantile, MultiStream, groupby
from .datatypes import SAMEntry, GFFEntry, Bed6, Interval, LocationEntry
from .io.strops import str_equal
from .util.cli import run_as_commandline
from .computation_graph import compute
from . import simulate
from . import arithmetics
from .genomic_data import Genome, GenomicArray, GenomicIntervals
from .plotting import plot
from .io.matrix_dump import Matrix
open = bnp_open


SAMBuffer = get_bufferclass_for_datatype(SAMEntry)
GFFBuffer = get_bufferclass_for_datatype(GFFEntry)
# Bed6Buffer = get_bufferclass_for_datatype(Bed6)

__all__ = ["EncodedArray", "EncodedRaggedArray",
           "KmerEncoder", "Minimizers", "PositionWeightMatrix", "mean",
           "bincount", "streamable", "histogram", "count_entries", "quantile",
           "BedBuffer", "VCFBuffer", "PhasedVCFMatrixBuffer", "VCFMatrixBuffer",
           "GfaSequenceBuffer",
           "TwoLineFastaBuffer", "FastQBuffer", "open_indexed", "groupby",
           "SAMBuffer", "GFFBuffer", "Bed6Buffer", "MultiLineFastaBuffer",
           "count_encoded", "DNAEncoding", "RNAENcoding", "AminoAcidEncoding"]


def set_backend(lib):
    def temp_insert(array, index, value):
        assert index == 0
        return lib.concatenate((lib.asarray([value], dtype=array.dtype), array))

    lib.insert = temp_insert

    nps.set_backend(lib)
    from npstructures.bitarray import BitArray
    from npstructures import RaggedArray
    from npstructures import RaggedShape, RaggedView

    if not hasattr(lib, "object_"):
        # hack for cupy
        lib.object_ = None

    from .encodings import set_backend as set_encoding_backend
    set_encoding_backend(lib)

    from .sequence import set_backend as set_sequence_backend
    set_sequence_backend(lib)

    from bionumpy.cupy_compatible.parser import CupyFileReader

    from .io import file_buffers
    file_buffers.np = lib

    from .io import parser
    parser.np = lib
    parser.NumpyFileReader = CupyFileReader

    from .io import files
    files.np = lib 
    files.NumpyFileReader = CupyFileReader

    from . import bnpdataclass
    bnpdataclass.np = lib

    from . import encoded_array
    encoded_array.np = lib
    encoded_array.get_NPSArray = lambda x: x

    from . import util
    util.np = lib

