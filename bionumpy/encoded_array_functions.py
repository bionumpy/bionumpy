from typing import List
import numpy as np
from npstructures import RaggedArray
from .encoded_array import EncodedArray, EncodedRaggedArray
from .encoded_array import EncodingException
from .encodings import BaseEncoding
from .encodings.identity_encoding import IdentityEncoding


def _list_of_encoded_arrays_as_encoded_ragged_array(array_list: List[EncodedArray]):
    assert all(isinstance(a, EncodedArray) for a in array_list)
    encoding = array_list[0].encoding
    assert all(a.encoding == encoding for a in array_list)
    data = np.concatenate([a.data for a in array_list])
    shape = [len(a) for a in array_list]
    return EncodedRaggedArray(EncodedArray(data, encoding), shape)


def _is_encoded(data):
    return isinstance(data, (EncodedArray, EncodedRaggedArray))


def as_encoded_array(s, target_encoding: "Encoding" = None) -> EncodedArray:
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
    if target_encoding.is_numeric() and type(s) in (np.ndarray, RaggedArray):
        return s
    # is already encoded if list and elements are encoded
    elif isinstance(s, list) and len(s) > 0 and isinstance(s[0], EncodedArray):
        return _list_of_encoded_arrays_as_encoded_ragged_array(s)

    return target_encoding.encode(s)


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
