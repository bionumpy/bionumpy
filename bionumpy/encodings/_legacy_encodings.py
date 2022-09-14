import numpy as np

class ThreeBitEncoding:
    reverse = np.array(
        [0, ord("A"), 0, ord("C"), ord("G"), 0, 0, ord("T")], dtype=np.uint8
    )
    complements = np.array([0, 7, 0, 4, 3, 0, 0, 1], dtype=np.uint8)

    @classmethod
    def encode(cls, sequence):
        return sequence & (2 ** 3 - 1)

    @classmethod
    def decode(cls, sequence):
        return cls.reverse[sequence]

    @classmethod
    def complement(cls, sequence):
        return cls.complements[sequence]


class BaseEncoding:
    """ Basic ACII byte encoding """

    complements = np.zeros(256, dtype=np.uint8)
    complements[[ord(c) for c in "ACGT"]] = [ord(c) for c in "TGCA"]
    complements[[ord(c) for c in "acgt"]] = [ord(c) for c in "tgca"]

    @classmethod
    def complement(cls, sequence):
        return cls.complements[sequence]

    @classmethod
    def from_string(cls, sequence):
        return np.array([ord(c) for c in sequence], dtype=np.uint8)

    @classmethod
    def encode(cls, sequence):
        """Identity"""
        return sequence

    @classmethod
    def decode(cls, sequence):
        """Identity"""
        return sequence

    @classmethod
    def to_string(cls, byte_sequence):
        return "".join(chr(b) for b in byte_sequence)


class ACTGTwoBitEncoding:
    letters = ["A", "C", "T", "G"]
    bitcodes = ["00", "01", "10", "11"]
    reverse = np.array([1, 3, 20, 7], dtype=np.uint8)
    _lookup_2bytes_to_4bits = np.zeros(256 * 256, dtype=np.uint8)
    _lookup_2bytes_to_4bits[
        256 * reverse[np.arange(4)[:, None]] + reverse[np.arange(4)]
    ] = np.arange(4)[:, None] * 4 + np.arange(4)
    _shift_4bits = 4 * np.arange(2, dtype=np.uint8)
    _shift_2bits = 2 * np.arange(4, dtype=np.uint8)

    @classmethod
    def convert_2bytes_to_4bits(cls, two_bytes):
        assert two_bytes.dtype == np.uint16, two_bytes.dtype
        return cls._lookup_2bytes_to_4bits[two_bytes]

    @classmethod
    def join_4bits_to_byte(cls, four_bits):
        return np.sum(four_bits << cls._shift_4bits, axis=1, dtype=np.uint8)

    @classmethod
    def complement(cls, char):
        complements = np.packbits([1, 0, 1, 0, 1, 0, 1, 0])
        dtype = char.dtype
        return (char.view(np.uint8) ^ complements).view(dtype)

    @classmethod
    def encode(cls, sequence):
        if sequence.size % 16 != 0:
            sequence = np.append(
                sequence, np.empty(16 - (sequence.size % 16), dtype=np.uint8)
            )

        assert sequence.dtype == np.uint8
        assert sequence.size % 4 == 0, sequence.size
        sequence = sequence & 31
        four_bits = cls.convert_2bytes_to_4bits(sequence.view(np.uint16))
        codes = cls.join_4bits_to_byte(four_bits.reshape(-1, 2))
        assert codes.dtype == np.uint8, codes.dtype
        return codes.flatten().view(np.uint8)

    @classmethod
    def from_string(cls, string):
        byte_repr = np.array([ord(c) for c in string], dtype=np.uint8)
        return cls.encode(byte_repr)

    @classmethod
    def to_string(cls, bits):
        byte_repr = cls.decode(bits)
        return "".join(chr(b) for b in byte_repr)

    @classmethod
    def decode(cls, sequence):
        assert sequence.dtype == np.uint8
        bit_mask = np.uint8(3)  # last two bits
        all_bytes = (sequence[:, None] >> cls._shift_2bits) & bit_mask
        return cls.reverse[all_bytes.flatten()] + 96


class ACTGEncoding:
    _lookup_byte_to_2bits = np.zeros(256, dtype=np.uint8)
    _lookup_byte_to_2bits[[97, 65]] = 0
    _lookup_byte_to_2bits[[99, 67]] = 1
    _lookup_byte_to_2bits[[116, 84]] = 2
    _lookup_byte_to_2bits[[103, 71]] = 3
    reverse = np.array([ord(c) for c in "ACTG"], dtype=np.uint8)

    @classmethod
    def encode(cls, bytes_array):
        return cls._lookup_byte_to_2bits[bytes_array]

    @classmethod
    def decode(cls, encoded):
        return cls.reverse[encoded]


class SimpleEncoding(ACTGTwoBitEncoding):
    _lookup_byte_to_2bits = np.zeros(256, dtype=np.uint8)
    _lookup_byte_to_2bits[[97, 65]] = 0
    _lookup_byte_to_2bits[[99, 67]] = 1
    _lookup_byte_to_2bits[[116, 84]] = 2
    _lookup_byte_to_2bits[[103, 71]] = 3

    _shift_2bits = 2 * np.arange(4, dtype=np.uint8)

    @classmethod
    def convert_byte_to_2bits(cls, one_byte):
        assert one_byte.dtype == np.uint8, one_byte.dtype
        return cls._lookup_byte_to_2bits[one_byte]

    @classmethod
    def join_2bits_to_byte(cls, two_bits_vector):
        return np.bitwise_or.reduce(two_bits_vector << cls._shift_2bits, axis=-1)

    @classmethod
    def encode(cls, sequence):
        assert sequence.dtype == np.uint8
        assert sequence.size % 4 == 0, sequence.size
        two_bits = cls.convert_byte_to_2bits(sequence)
        codes = cls.join_2bits_to_byte(two_bits.reshape(-1, 4))
        return codes.flatten()


def twobit_swap(number):
    dtype = number.dtype
    byte_lookup = np.zeros(256, dtype=np.uint8)
    power_array = 4 ** np.arange(4)
    rev_power_array = power_array[::-1]
    for two_bit_string in product([0, 1, 2, 3], repeat=4):
        byte_lookup[np.sum(power_array * two_bit_string)] = np.sum(
            rev_power_array * two_bit_string
        )
    new_bytes = byte_lookup[number.view(np.uint8)]
    return new_bytes.view(dtype).byteswap()


def alphabet_encoding(_alphabet, name):
    _alphabet = np.asanyarray(_alphabet)
    _lookup = np.zeros(np.max(_alphabet) + 1, dtype=np.uint8)
    _lookup[_alphabet] = np.arange(_alphabet.size)

    class cls:
        alphabet = _alphabet
        lookup = _lookup

        @classmethod
        def encode(cls, byte_array):
            return cls.lookup[byte_array]

        @classmethod
        def decode(cls, encoded):
            return cls.alphabet[encoded]

    cls.__name__ = name
    return cls


AminoAcidEncoding = alphabet_encoding(
    [65, 67, 68, 69, 70, 71, 72, 73, 75, 76, 77, 78, 80, 81, 82, 83, 84, 86, 87, 89],
    "AminoAcidEncoding",
)
