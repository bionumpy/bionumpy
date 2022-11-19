from abc import abstractmethod
import numpy as np


class Encoding:
    @abstractmethod
    def encode(self, *args, **kwargs):
        return NotImplemented

    @abstractmethod
    def get_labels(self):
        pass


class OneToOneEncoding(Encoding):
    @abstractmethod
    def encode(self, ascii_codes):
        return NotImplemented

    @abstractmethod
    def decode(self, encoded):
        return NotImplemented


class ASCIIEncoding(OneToOneEncoding):
    def encode(self, ascii_codes):
        return ascii_codes

    def decode(self, encoded):
        return encoded

    def __repr__(self):
        return "ASCIIEncoding()"

    def __hash__(self):
        return hash(repr(self))


class NumericEncoding(OneToOneEncoding):
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
