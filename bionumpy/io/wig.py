from .delimited_buffers import DelimitedBuffer, DelimitedBufferWithInernalComments
from ..datatypes import BedGraph
from .strops import split
from .exceptions import FormatException
from ..encodings.alphabet_encoding import DigitEncoding
from ..encodings.exceptions import EncodingError
from ..encoded_array import BaseEncoding, EncodedArray, as_encoded_array
import numpy as np


class WigBuffer(DelimitedBufferWithInernalComments):
    dataclass = BedGraph
    DELIMITER = '\t'

