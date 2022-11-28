from typing import List
import numpy as np
from npstructures import RaggedArray
from npstructures.testing import assert_raggedarray_equal

from .encoded_array import EncodedArray, EncodedRaggedArray
from .encoded_array import EncodingException, IncompatibleEncodingsException
from .encodings import Encoding, BaseEncoding
from .encodings.identity_encoding import IdentityEncoding


def encode_string(s: str, target_encoding):
    s = EncodedArray([ord(c) for c in s], IdentityEncoding)
    s = _encode_base_encoded_array(s, target_encoding)
    return s


def list_of_encoded_arrays_as_encoded_ragged_array(array_list: List[EncodedArray]):
    assert all(isinstance(a, EncodedArray) for a in array_list)
    encoding = array_list[0].encoding
    assert all(a.encoding == encoding for a in array_list)
    data = np.concatenate([a.data for a in array_list])
    shape = [len(a) for a in array_list]
    return EncodedRaggedArray(EncodedArray(data, encoding), shape)


def encode_list_of_strings(s: str, target_encoding):
    s = EncodedRaggedArray(
        EncodedArray([ord(c) for ss in s for c in ss], IdentityEncoding),
        [len(ss) for ss in s])
    return ragged_array_as_encoded_array(s, target_encoding)


def _is_encoded(data):
    return isinstance(data, (EncodedArray, EncodedRaggedArray))


def as_encoded_array(s, target_encoding: Encoding = None) -> EncodedArray:
    """Main function used to create encoded arrays from e.g. strings orl lists.
    Can be called with already encoded arrays, and will then do nothing.

    If input is `str` or `List[str]` objects, creates an `EncodedArray` or `EncodedRaggedArray`
    object from them with the given encoding.

    If the input is an `EncodedArray` or `EncodedRaggedArray` AND input is BaseEncoded,
    encode the input to the `target_encoding` if possible. If `target_encoding` is None, nothing is done.

    Raw encoded data as `np.ndarray` objects should not be passed to this function. If you have
    already encoded data in `np.ndarray` objects, use the `EncodedArray.__init__`directly

    Parameters
    ----------
    s : str/List[str]/EnocedArray/EncodedRaggedArray
        The data to be represented in an EncodedArray
    target_encoding : Encoding
        The encoding to use in the resulting EncodedArray

    Returns
    -------
    EncodedArray
        Encoded data in an EncodedArray

    default target encoding None:
    if None: encode as base encoding if it is not encoded
    if already encoded: do nothing
    this function is not for changing encoding on stuff

    """
    if isinstance(s, (EncodedArray, EncodedRaggedArray)):
        if target_encoding is None or s.encoding == target_encoding:
            return s
        else:
            if not s.encoding.is_base_encoding():
                raise EncodingException("Trying to encode already encoded array with encoding %s to encoding %s. "
                                        "This is not supported. Use the change_encoding function." % (
                    s.encoding, target_encoding))
    elif target_encoding is None:
        target_encoding = BaseEncoding

    # if numeric encoding and already np-array, this is already encoded
    if target_encoding.is_numeric() and isinstance(s, (np.ndarray, RaggedArray)):
        return s

    if isinstance(s, list) and len(s) > 0 and isinstance(s[0], EncodedArray):
        return list_of_encoded_arrays_as_encoded_ragged_array(s)
    elif isinstance(s, (RaggedArray, EncodedRaggedArray)):
        #assert hasattr(target_encoding, "_encode"), target_encoding
        new = target_encoding.encode(s)
        #assert isinstance(new, RaggedArray), "%s, %s" % (type(new), type(s))
        #old = ragged_array_as_encoded_array(s, target_encoding)
        #assert np.all(old.ravel() == new.ravel()), "%s, %s != %s, %s, %s, %s, %s" % (target_encoding, old, new, type(old), type(new), repr(old), repr(new))
        #assert old.shape == new.shape, "%s != %s, %s, %s, %s, %s" % (old, new, type(old), type(new), repr(old), repr(new))
        #assert_raggedarray_equal(old, new), "%s != %s, %s, %s, %s, %s" % (old, new, type(old), type(new), repr(old), repr(new))
        return new
    elif isinstance(s, np.ndarray):
        return np_array_as_encoded_array(s, target_encoding)
    else:
        if _is_encoded(s) and target_encoding == s.encoding:
            return s
        else:
            return target_encoding.encode(s)


def _encode_encoded_array(encoded_array, target_encoding):
    """Encode an encoded array with a new encoding.
    Encoded array should either be BaseEncoded or have target_encoding"""
    assert isinstance(encoded_array, EncodedArray), (encoded_array, repr(encoded_array), type(encoded_array))
    if encoded_array.encoding == target_encoding:
        return encoded_array

    if encoded_array.encoding.is_base_encoding():
        encoded_array = _encode_base_encoded_array(encoded_array, target_encoding)
    elif target_encoding.is_base_encoding():
        encoded_array = EncodedArray(encoded_array.encoding.decode(encoded_array.data), BaseEncoding)
    else:
        raise IncompatibleEncodingsException("Can only encode EncodedArray with BaseEncoding or target encoding. "
                                             "Base encoding is %s, target encoding is %s"
                                            % (str(encoded_array.encoding), str(target_encoding)))
    return encoded_array


def _encode_base_encoded_array(encoded_array, target_encoding):
    assert encoded_array.encoding.is_base_encoding()
    encoded_array = target_encoding.encode(encoded_array.data)
    if target_encoding.is_numeric():
        encoded_array = encoded_array
    else:
        encoded_array = EncodedArray(encoded_array, target_encoding)
    return encoded_array


def np_array_as_encoded_array(s, target_encoding):
    assert hasattr(target_encoding, "is_numeric")
    return target_encoding.encode(s)


def ragged_array_as_encoded_array(s, target_encoding):
    data = as_encoded_array(s.ravel(), target_encoding)
    if isinstance(data, EncodedArray):
        return EncodedRaggedArray(data, s.shape)
    return RaggedArray(data, s.shape)


def from_encoded_array(encoded_array: EncodedArray) -> str:
    """Convert data in an `EncodedArray`/`EncodedRaggedArray into `str`/`List[str]`

    Unlike the `EncodedArray.__str__` this will convert all the data into strings

    Parameters
    ----------
    encoded_array : EncodedArray

    Returns
    -------
    str
        Full string representation

    Examples
    --------
    5

    """
    if isinstance(encoded_array, EncodedRaggedArray):
        return [from_encoded_array(row) for row in encoded_array]
    else:
        return "".join(chr(c) for c in encoded_array.encoding.decode(encoded_array).raw())


def change_encoding(encoded_array, new_encoding):
    assert isinstance(encoded_array, (EncodedArray, EncodedRaggedArray)), \
        "Can only change encoding of EncodedArray or EncodedRaggedArray"

    new_data = new_encoding.encode(
        encoded_array.encoding.decode(encoded_array.ravel())
    )

    if isinstance(encoded_array, EncodedArray):
        return EncodedArray(new_data, new_encoding)
    elif isinstance(encoded_array, EncodedRaggedArray):
        return EncodedRaggedArray(EncodedArray(new_data, new_encoding), encoded_array.shape)
