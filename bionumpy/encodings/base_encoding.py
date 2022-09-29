from abc import abstractmethod
import numpy as np


class Encoding:

    @abstractmethod
    def encode(self, ascii_codes):
        return NotImplemented

    @abstractmethod
    def decode(self, encoded):
        return NotImplemented


class ASCIIEncoding(Encoding):
    def encode(self, ascii_codes):
        return ascii_codes

    def decode(self, encoded):
        return encoded


class NumericEncoding(Encoding):
    pass


class CigarEncoding(NumericEncoding):
    @classmethod
    def encode(cls, data):
        return np.asarray(data)

    @classmethod
    def decode(cls, data):
        return np.asarray(data)


BaseEncoding = ASCIIEncoding()
