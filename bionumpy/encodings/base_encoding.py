from abc import abstractmethod
import numpy as np
from npstructures import RaggedArray

from bionumpy.encodings.identity_encoding import IdentityEncoding


class Encoding:
    @abstractmethod
    def encode(self, *args, **kwargs):
        return NotImplemented

    @abstractmethod
    def get_labels(self):
        pass

    def is_base_encoding(self):
        return False

    def is_one_to_one_encoding(self):
        return False


class OneToOneEncoding(Encoding):
    def encode(self, data):
        from ..encoded_array import EncodedArray, EncodedRaggedArray
        from ..encoded_array_functions import list_of_encoded_arrays_as_encoded_ragged_array
        assert hasattr(self, "_encode"), "Missing implementation of _encode for %s" % self
        if isinstance(data, (EncodedArray, EncodedRaggedArray)):
            assert data.encoding.is_base_encoding(), "Data is already encoded. " \
                                                     "Can only encode already encoded data if it is base encoded."
            data = data.raw()  # input from here is always "raw"

        if isinstance(data, str):
            out = self._encode_string(data)
        elif isinstance(data, list) and len(data) > 0 and isinstance(data[0], EncodedArray):
            out = list_of_encoded_arrays_as_encoded_ragged_array(data)
        elif isinstance(data, list):
            out = self._encode_list_of_strings(data)
        elif isinstance(data, RaggedArray):
            print("IS RAGGED ARRAY")
            r = self._ragged_array_as_encoded_array(data)
            assert isinstance(r, (EncodedRaggedArray, RaggedArray))
            return r
        elif isinstance(data, np.ndarray):
            if isinstance(self, NumericEncoding):
                out = self._encode(data)
            else:
                out = EncodedArray(self._encode(data), self)
        else:
            assert False

        return out

    def _encode_list_of_strings(self, s: str):
        from ..encoded_array import EncodedArray, EncodedRaggedArray
        s = EncodedRaggedArray(
            EncodedArray([ord(c) for ss in s for c in ss], IdentityEncoding()),
            [len(ss) for ss in s])
        return self._ragged_array_as_encoded_array(s)

    def _ragged_array_as_encoded_array(self, s):
        from ..encoded_array import EncodedArray, EncodedRaggedArray

        data = self.encode(s.ravel())
        if isinstance(data, EncodedArray):
            return EncodedRaggedArray(data, s.shape)

        print("_ragged array as encoded array")
        print(s)

        return RaggedArray(data, s.shape)

        """
        if isinstance(s, EncodedRaggedArray):
            s = s.raw()

        if isinstance(self, NumericEncoding):
            data = self._encode(s.ravel())
            out_class = RaggedArray
        else:
            out_class = EncodedRaggedArray
            data = EncodedArray(self._encode(s.ravel()), self)

        return out_class(data, s.shape)
        """

    def _encode_string(self, string: str):
        from ..encoded_array import EncodedArray, EncodedRaggedArray
        s = EncodedArray([ord(c) for c in string], IdentityEncoding())
        s = self._encode_base_encoded_array(s)
        return s

    def _encode_base_encoded_array(self, encoded_array):
        from ..encoded_array import EncodedArray, EncodedRaggedArray
        assert encoded_array.encoding.is_base_encoding()
        encoded_array = self._encode(encoded_array.data)
        if hasattr(self, "is_numeric"):
            encoded_array = encoded_array
        else:
            encoded_array = EncodedArray(encoded_array, self)
        return encoded_array

    def decode(self, data):
        from ..encoded_array import EncodedArray, EncodedRaggedArray
        if not hasattr(self, "_decode"):
            raise Exception("Missing implementation of _decode for %s" % self)

        if isinstance(data, int):
            return EncodedArray(self._decode(np.atleast_1d(data)), self)
        elif isinstance(data, np.ndarray):
            assert isinstance(self, NumericEncoding), "%s" % data
            return self._decode(data)
        elif isinstance(data, EncodedRaggedArray):
            return EncodedRaggedArray(
                EncodedArray(self._decode(data.raw().ravel()), BaseEncoding), data.shape)
        elif isinstance(data, RaggedArray):
            assert isinstance(self, NumericEncoding), "%s" % data
            return RaggedArray(self._decode(data.ravel()), data.shape)
        elif isinstance(data, EncodedArray):
            return EncodedArray(self._decode(data.raw()), BaseEncoding)
        else:
            raise Exception("Not able to decode %s with %s" % (data, self))


    def is_one_to_one_encoding(self):
        return True


class ASCIIEncoding(OneToOneEncoding):
    def _encode(self, ascii_codes):
        return ascii_codes

    def _decode(self, encoded):
        return encoded

    def __repr__(self):
        return "ASCIIEncoding()"

    def __hash__(self):
        return hash(repr(self))

    def is_base_encoding(self):
        return True


class NumericEncoding(OneToOneEncoding):
    is_numeric = True


"""
class CigarEncoding(NumericEncoding):
    @classmethod
    def encode(cls, data):
        return np.asarray(data)

    @classmethod
    def decode(cls, data):
        return np.asarray(data)
"""


BaseEncoding = ASCIIEncoding()


def get_base_encodings():
    return [BaseEncoding]  # TODO: add other encodings
