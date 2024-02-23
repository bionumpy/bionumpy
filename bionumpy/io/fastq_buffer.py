import dataclasses
from typing import List

import numpy as np

from .dump_csv import get_column
from ..encoded_array import EncodedRaggedArray, as_encoded_array, EncodedArray, BaseEncoding, change_encoding
from ..datatypes import SequenceEntryWithQuality
from ..encodings import QualityEncoding
from ..io.exceptions import FormatException
from .one_line_buffer import OneLineBuffer


class FastQBuffer(OneLineBuffer):
    HEADER = "@"
    n_lines_per_entry = 4
    dataclass = SequenceEntryWithQuality
    _line_offsets = (1, 0, 0, 0)
    _empty_lines = [2]

    def get_text_field_by_number(self, i: int) -> EncodedRaggedArray:
        if i == 2:
            return self._buffer_extractor.get_field_by_number(3)
        return super().get_text_field_by_number(i)

    def get_field_by_number(self, i: int, t: type=None):
        if i == 2:
            return QualityEncoding.encode(self.get_text_field_by_number(i))
        else:
            return super().get_field_by_number(i, t)

    def get_data(self):
        seq_entry = super().get_data()
        quality = self.get_field_by_number(2, QualityEncoding)
        return SequenceEntryWithQuality(
            seq_entry.name, seq_entry.sequence,quality)

    @classmethod
    def _validate(cls, data, new_lines):
        super()._validate(data, new_lines)
        n_lines_per_entry = cls.n_lines_per_entry
        if np.any(data[new_lines[1::n_lines_per_entry] + 1] != "+"):
            entry_number = np.flatnonzero(data[new_lines[1::n_lines_per_entry] + 1] != "+")[0]
            line_number = 2 + entry_number * n_lines_per_entry
            raise FormatException(f"Expected '+' at third line of entry in {data}", line_number=line_number)

    @classmethod
    def join_fields(cls, fields: List[EncodedRaggedArray]):
        plus_line = as_encoded_array(['+'] * len(fields[0]))
        return super().join_fields(fields[:2]+[plus_line]+fields[2:])

    @classmethod
    def from_data(cls, entries):
        name_field = get_column(entries.name, dataclasses.fields(entries)[0].type)
        quality_field = EncodedRaggedArray(EncodedArray(QualityEncoding.decode(entries.quality.ravel()), BaseEncoding),
                                           entries.quality.shape)
        sequence_field = change_encoding(entries.sequence, BaseEncoding) if entries.sequence.encoding != BaseEncoding else entries.sequence


        return cls.join_fields([name_field, sequence_field,
                                quality_field])
