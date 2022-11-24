from abc import abstractmethod
import numpy as np


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
    @abstractmethod
    def encode(self, ascii_codes):
        return NotImplemented

    def encode(self, data):
        from ..encoded_array import EncodedArray, EncodedRaggedArray
        assert hasattr(self, "_encode"), "Missing implementation of _encode for %s" % self


        if isinstance(data, (EncodedArray, EncodedRaggedArray)):
            assert data.encoding.is_base_encoding(), "Data is already encoded. " \
                                                     "Can only encode already encoded data if it is base encoded."

        return self._encode(data)

    def decode(self, data):
        if hasattr(self, "_decode"):
            return self._decode(data)

        raise Exception("Missing implementation of _decode for %s" % self)

    def is_one_to_one_encoding(self):
        return True


class ASCIIEncoding(OneToOneEncoding):
    def encode(self, ascii_codes):
        return ascii_codes

    def decode(self, encoded):
        return encoded

    def __repr__(self):
        return "ASCIIEncoding()"

    def __hash__(self):
        return hash(repr(self))

    def is_base_encoding(self):
        return True


class NumericEncoding(OneToOneEncoding):
    is_numeric = True
    pass


class CigarEncoding(NumericEncoding):
    @classmethod
    def encode(cls, data):
        return np.asarray(data)

    @classmethod
    def decode(cls, data):
        return np.asarray(data)


BaseEncoding = ASCIIEncoding()


def get_base_encodings():
    return [BaseEncoding]  # TODO: add other encodings
