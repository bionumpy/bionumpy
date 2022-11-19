"""Top-level package for bionumpy."""

__author__ = """Knut Rand"""
__email__ = "knutdrand@gmail.com"
__version__ = '0.2.5'

import npstructures as nps

from .io.files import bnp_open as open
from .io.files import count_entries, open_indexed
from .io.multiline_buffer import MultiLineFastaBuffer
from .io.file_buffers import (TwoLineFastaBuffer, FastQBuffer)
from .io.delimited_buffers import (BedBuffer, VCFBuffer, VCFMatrixBuffer,
                                GfaSequenceBuffer, get_bufferclass_for_datatype)
from .encodings.alphabet_encoding import (DNAEncoding, RNAENcoding, AminoAcidEncoding)
from .encoded_array import EncodedArray, EncodedRaggedArray, as_encoded_array
from .kmers import KmerEncoding
from .minimizers import Minimizers
from .position_weight_matrix import PositionWeightMatrix
from .counter import count_encoded
from .streams import mean, bincount, histogram, streamable, quantile, MultiStream, groupby
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
    def temp_insert(array, index, value):
        assert index == 0
        return lib.concatenate((lib.asarray([value], dtype=array.dtype), array))

    lib.insert = temp_insert

    import sys

    nps.set_backend(lib)

    from .encodings import set_backend as set_encoding_backend
    set_encoding_backend(lib)
    
    from npstructures.bitarray import BitArray
    from npstructures import RaggedArray
    from npstructures import RaggedShape, RaggedView

    #from bionumpy.cupy_compatible.sequences import CPSequence, CPSequences
    #from bionumpy.cupy_compatible.parser import CpBufferStream
    #from bionumpy.cupy_compatible.file_buffers import CPTwoLineFastaBuffer
    #from bionumpy.cupy_compatible.encodings.alphabet_encoding import CPAlphabetEncoding
    #from bionumpy.cupy_compatible.encodings.alphabet_encoding import CPACTGEncoding
    #from bionumpy.cupy_compatible.encodings.alphabet_encoding import CPACTGnEncoding
    #from bionumpy.cupy_compatible.encodings.alphabet_encoding import CPAminoAcidEncoding

    from bionumpy.cupy_compatible.parser import CupyFileReader

    #sys.modules[__name__].Sequence = CPSequence
    #sys.modules[__name__].Sequences = CPSequences

    from .io import file_buffers
    file_buffers.np = lib
    #file_buffers.TwoLineFastaBuffer = CPTwoLineFastaBuffer
    file_buffers.RaggedArray = RaggedArray
    file_buffers.RaggedShape = RaggedShape
    file_buffers.RaggedView = RaggedView

    from .io import parser
    parser.np = lib
    parser.NumpyFileReader = CupyFileReader

    from .io import files
    files.np = lib 
    files.NumpyFileReader = CupyFileReader
    #files.NpBufferStream = CpBufferStream
    #files.buffer_types[".fa"] = CPTwoLineFastaBuffer
    #files.buffer_types[".fasta"] = CPTwoLineFastaBuffer

    from . import bnpdataclass
    bnpdataclass.np = lib
    bnpdataclass.RaggedArray = RaggedArray

    from . import encoded_array
    encoded_array.np = lib
    encoded_array.RaggedArray = RaggedArray
    encoded_array.get_NPSArray = lambda x: x


    from . import kmers
    #kmers.ACTGEncoding = ACTGEncoding
    kmers.BitArray = BitArray

    from . import util
    util.np = lib
    util.RaggedArray = RaggedArray
