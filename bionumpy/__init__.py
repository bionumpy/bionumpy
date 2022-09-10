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
