import re
from ..encoded_array import as_encoded_array
from ..encoded_array import EncodedArray, EncodedRaggedArray
from ..encodings import BaseEncoding
import numpy as np

from ..string_array import as_string_array


def match_regexp(encoded_array, regexp: str):
    encoded_array = as_encoded_array(encoded_array)
    s = re.findall(regexp, encoded_array.to_string())
    return EncodedRaggedArray(
        EncodedArray(np.frombuffer(bytes(''.join(s), 'ascii'), dtype=np.uint8), BaseEncoding),
        [len(ss) for ss in s])

def match_regexp_string_array(encoded_array, regexp: str):
    encoded_array = as_encoded_array(encoded_array)
    s = re.findall(regexp, encoded_array.to_string())
    return as_string_array(s)

