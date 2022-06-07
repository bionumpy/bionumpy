from abc import abstractmethod


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


BaseEncoding = ASCIIEncoding()
