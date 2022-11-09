"""Top-level package for bionumpy."""

__author__ = """Knut Rand"""
__email__ = "knutdrand@gmail.com"
__version__ = '0.2.3'

import npstructures as nps

from .io.files import bnp_open as open
from .io.files import count_entries, open_indexed
from .io.multiline_buffer import MultiLineFastaBuffer
from .io.file_buffers import (TwoLineFastaBuffer, FastQBuffer)
from .io.delimited_buffers import (BedBuffer, VCFBuffer, VCFMatrixBuffer,
                                GfaSequenceBuffer, get_bufferclass_for_datatype)
from .encodings.alphabet_encoding import (DNAEncoding, RNAENcoding, AminoAcidEncoding)
from .encoded_array import EncodedArray, EncodedRaggedArray, as_encoded_array
from .groupby import groupby
from .kmers import KmerEncoding
from .minimizers import Minimizers
from .position_weight_matrix import PositionWeightMatrix
from .counter import count_encoded
from .streams import mean, bincount, histogram, streamable, quantile, MultiStream
from .datatypes import SAMEntry, GFFEntry, Bed6


SAMBuffer = get_bufferclass_for_datatype(SAMEntry)
GFFBuffer = get_bufferclass_for_datatype(GFFEntry)
Bed6Buffer = get_bufferclass_for_datatype(Bed6)

__all__ = ["EncodedArray", "as_encoded_array", "EncodedRaggedArray",
           "KmerEncoding", "Minimizers", "PositionWeightMatrix", "mean",
           "bincount", "streamable", "histogram", "count_entries", "quantile",
           "BedBuffer", "VCFBuffer", "VCFMatrixBuffer", "GfaSequenceBuffer",
           "TwoLineFastaBuffer", "FastQBuffer", "open_indexed", "groupby",
           "SAMBuffer", "GFFBuffer", "Bed6Buffer", "MultiLineFastaBuffer",
           "count_encoded", "DNAEncoding", "RNAENcoding", "AminoAcidEncoding"]


def set_backend(lib):
    import sys

    nps.set_backend(lib)
    
    from npstructures.bitarray import BitArray
    from npstructures import RaggedArray
    from npstructures import RaggedShape, RaggedView

    from bionumpy.cupy_compatible.sequences import CPSequence, CPSequences
    from bionumpy.cupy_compatible.parser import CpBufferStream
    from bionumpy.cupy_compatible.file_buffers import CPTwoLineFastaBuffer
    #from bionumpy.cupy_compatible.encodings.alphabet_encoding import CPAlphabetEncoding
    #from bionumpy.cupy_compatible.encodings.alphabet_encoding import CPACTGEncoding
    #from bionumpy.cupy_compatible.encodings.alphabet_encoding import CPACTGnEncoding
    #from bionumpy.cupy_compatible.encodings.alphabet_encoding import CPAminoAcidEncoding

    sys.modules[__name__].Sequence = CPSequence
    sys.modules[__name__].Sequences = CPSequences

    from . import file_buffers
    file_buffers.np = lib
    file_buffers.TwoLineFastaBuffer = CPTwoLineFastaBuffer
    #file_buffers.RaggedShape = RaggedShape
    #file_buffers.RaggedView = RaggedView

    from . import parser
    parser.np = lib

    from . import files
    files.np = lib 
    files.NpBufferStream = CpBufferStream
    files.buffer_types[".fa"] = CPTwoLineFastaBuffer
    files.buffer_types[".fasta"] = CPTwoLineFastaBuffer

    from .encodings import set_backend as set_encoding_backend
    set_encoding_backend(lib)

    from . import kmers
    #kmers.ACTGEncoding = ACTGEncoding
    kmers.BitArray = BitArray

    from . import util
    util.RaggedArray = RaggedArray
