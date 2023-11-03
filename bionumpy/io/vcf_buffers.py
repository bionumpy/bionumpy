import dataclasses
import logging
from typing import List, Optional

import numpy as np
from npstructures import RaggedView, RaggedArray
from npstructures.raggedshape import RaggedView2

from .exceptions import FormatException
from .file_buffers import TextBufferExtractor, FileBuffer
from .vcf_header import parse_header
from ..bnpdataclass.bnpdataclass import narrow_type
from ..bnpdataclass.lazybnpdataclass import create_lazy_class, ItemGetter
from ..encoded_array import EncodedArray, EncodedRaggedArray, as_encoded_array
from ..bnpdataclass import BNPDataClass, make_dataclass
from ..datatypes import VCFEntry, VCFGenotypeEntry, PhasedVCFGenotypeEntry, PhasedVCFHaplotypeEntry
from ..encodings.vcf_encoding import GenotypeRowEncoding, PhasedGenotypeRowEncoding, PhasedHaplotypeRowEncoding
from .delimited_buffers import DelimitedBuffer

logger = logging.getLogger(__name__)


class GenotypeBufferExtractor(TextBufferExtractor):
    def get_genotypes(self):
        indices = self._starts[..., None] + np.arange(3)
        return self._data[indices]


class NamedBufferExtractor(TextBufferExtractor):
    def __init__(self, data, field_starts, field_lens, names):
        super().__init__(data, field_starts, field_lens=field_lens)
        self._names = names

    @classmethod
    def concatenate(cls, buffers: List['NamedBufferExtractor']):
        sizes = np.array([b._data.size for b in buffers])
        offsets = np.insert(np.cumsum(sizes), 0, 0)
        data = np.concatenate([b._data for b in buffers])
        starts = np.concatenate([b._field_starts + offset for b, offset in zip(buffers, offsets)])
        lens = np.concatenate([b._field_lens for b in buffers])
        return cls(data, starts, field_lens=lens, names=buffers[0]._names)

    def __getitem__(self, idx):
        return self.__class__(self._data, field_starts=self._field_starts[idx], field_lens=self._field_lens[idx], names=self._names)

    def get_field_by_number(self, number, keep_sep=False):
        name = self._names[number]
        return self.get_field_by_name(name, keep_sep=keep_sep)

    def has_field_number(self, number):
        name = self._names[number]
        return self.has_field_name(name)

    def has_field_name(self, name):
        starts = self._field_starts.ravel()
        mask = self._field_lens.ravel() == len(name)
        if np.any(mask):
            array = RaggedArray(self._data,
                                RaggedView2(starts[mask], np.full(mask.sum(), len(name)))).to_numpy_array()
            mask[mask] = (array == name).all(axis=-1)
        return RaggedArray(mask, self._field_starts.shape).any(axis=1)

    def get_field_by_name(self, name, keep_sep=False):
        assert name in self._names, (name, self._names)
        mask = self.has_field_mask(name)
        n_entries = len(self._field_starts)
        if not np.any(mask):
            logger.warning(f"Field: {name} not found in buffer")
            if keep_sep:
                return EncodedRaggedArray(as_encoded_array(';'*n_entries), np.ones(n_entries, dtype=int))
            return EncodedRaggedArray(as_encoded_array(''), np.zeros(n_entries, dtype=int))
        reshaped_mask = RaggedArray(mask, self._field_starts.shape)
        line_sums = reshaped_mask.sum(axis=-1)
        if np.any(line_sums > 1):
            raise FormatException(f"Field: {name} found multiple times in buffer", line_number=np.flatnonzero(line_sums > 1)[0])
        present_mask = reshaped_mask.any(axis=-1)#line_sums.astype(bool)#

        field_starts = self._field_starts.ravel()[mask] + len(name) + 1
        lens = self._field_lens.ravel()[mask]-len(name)-1# - field_starts
        if keep_sep:
            lens += 1
        starts = np.zeros(n_entries, dtype=int)
        starts[present_mask] = field_starts
        starts = np.maximum.accumulate(starts)
        #starts = np.maximum.accumulate(np.where(present_mask, field_starts, 0))
        all_lens = np.zeros(n_entries, dtype=int)
        all_lens[present_mask] = lens
        text = EncodedRaggedArray(self._data, RaggedView2(starts, all_lens))
        return text

    def has_field_mask(self, name):
        line_len = len(name) + 1
        starts = self._field_starts.ravel()
        n_ignored_fields = 0
        while starts[len(starts) - n_ignored_fields - 1] + line_len >= self._data.size:
            n_ignored_fields += 1
        starts = starts[:len(starts) - n_ignored_fields]
        e = EncodedRaggedArray(
            self._data,
            RaggedView2(starts, [line_len] * len(starts)))
        flat_e = e.ravel()
        mask = (flat_e.reshape(-1, line_len) == name + '=').all(axis=-1)
        if n_ignored_fields:
            mask = np.append(mask, np.full(n_ignored_fields, False))
        return mask


class InfoBuffer(DelimitedBuffer):
    def __init__(self, buffer_extractor: NamedBufferExtractor, dataclass: BNPDataClass):
        self._buffer_extractor = buffer_extractor
        self._dataclass = dataclass
        self._is_validated = True

    @classmethod
    def concatenate(cls, buffers: List['InfoBuffer']):
        extractor = buffers[0]._buffer_extractor.concatenate([b._buffer_extractor for b in buffers])
        return cls(extractor, buffers[0]._dataclass)

    def __getitem__(self, idx):
        return self.__class__(self._buffer_extractor[idx], self._dataclass)

    def _get_field_by_number(self, col_number, field_type):
        if field_type == bool:
            return self._buffer_extractor.has_field_number(col_number)
        return super()._get_field_by_number(col_number, field_type)

    def _validate(self):
        pass


def create_info_dataclass(header_data):
    if not header_data:
        return str
    header = parse_header(header_data)
    is_list = lambda val: (val['Number'] is None) or (val['Number'] > 1)
    is_int_list = lambda val: (val['Type'] == Optional[int]) and is_list(val)
    convert_type = lambda val: List[int] if is_int_list(val) else (str if is_list(val) else val['Type'])
    info_fields = [(key, convert_type(val)) for key, val in header.INFO.items()]
    dc = make_dataclass(info_fields, "InfoDataclass")
    return dc


class VCFBuffer(DelimitedBuffer):
    dataclass = VCFEntry
    lazy_dataclass = create_lazy_class(dataclass)
    _info_dataclass = None
    _vcf_data_class = None
    info_cache = {}
    vcfentry_cache = {}


    @property
    def actual_dataclass(self):
        return self.vcf_data_class

    def _get_field_by_number(self, field_nr: int, field_type: type = object):
        if field_nr == 7:
            return self._get_info_field()
        val = super()._get_field_by_number(field_nr, field_type)
        if field_nr == 1:
            val -= 1
        return val

    @classmethod
    def from_data(cls, data: BNPDataClass) -> "DelimitedBuffer":
        data = dataclasses.replace(data, position=data.position + 1)
        return super().from_data(data)

    @property
    def info_dataclass(self):
        if self._info_dataclass is None:
            self._info_dataclass = self._make_info_dataclass()
        return self._info_dataclass

    @property
    def vcf_data_class(self):
        if self._vcf_data_class is None:
            self._vcf_data_class = self._make_vcf_dataclass()
        return self._vcf_data_class

    def _get_info_field(self):
        text = self._buffer_extractor.get_field_by_number(7)
        if (not self._header_data) or ('##INFO' not in self._header_data):
            logger.warning('No header data found. Cannot parse INFO field. Returning as string')
            return text
        delimiters = np.flatnonzero(text.ravel() == ';') + 1
        offsets = np.insert(np.cumsum(text.lengths), 0, 0)
        all_delimiters = np.sort(np.concatenate([delimiters, offsets]), kind='mergesort')
        delimiter_offsets = np.searchsorted(all_delimiters, offsets)
        dl_lens = np.diff(delimiter_offsets)
        starts = RaggedArray(all_delimiters[:-1].copy(), dl_lens)
        ends = RaggedArray(all_delimiters[1:], dl_lens)
        # starts[:, :-1] = starts[:, :-1] + 1
        ends[:, :-1] = ends[:, :-1] - 1
        dataclass = self.info_dataclass
        lens = ends - starts
        extractor = NamedBufferExtractor(
            text.ravel(),
            starts,
            lens,
            [f.name for f in dataclasses.fields(dataclass)])
        buf = InfoBuffer(extractor, dataclass)
        item_getter = ItemGetter(buf, dataclass)
        return self._lazy_info_class(item_getter)
        return create_lazy_class(dataclass)(item_getter)

    @property
    def _lazy_info_class(self):
        return self.__class__.info_cache[self.header_data][1]

    @property
    def _lazy_vcf_class(self):
        if self.header_data not in self.__class__.vcfentry_cache:
            self._make_vcf_dataclass()
        return self.__class__.vcfentry_cache[self.header_data][1]

    def _make_vcf_dataclass(self):
        vcfentry_cache = self.__class__.vcfentry_cache
        dataclass = self.dataclass
        header_data = self.header_data
        if header_data in vcfentry_cache:
            return vcfentry_cache[header_data][0]
        info_class = str if not header_data else self.info_dataclass
        vcf_entry = narrow_type(dataclass, 'info', info_class)
        lc = create_lazy_class(vcf_entry)
        vcfentry_cache[header_data] = (vcf_entry, lc)
        return vcfentry_cache[header_data][0]

    @classmethod
    def modify_class_with_header_data(cls, header_data):
        if not header_data or '##INFO' not in header_data:
            return cls
        info_class = create_info_dataclass(header_data)
        new_dataclass = narrow_type(cls.dataclass, 'info', info_class)
        new_lazy_class = create_lazy_class(new_dataclass)

        class ModifiedClass(cls):
            _header_data = header_data
            dataclass = new_dataclass
            lazy_class = new_lazy_class

        ModifiedClass.__name__ = cls.__name__+'H'
        ModifiedClass.__qualname__ = cls.__qualname__+'H'

        return ModifiedClass

    def _make_info_dataclass(self):
        info_cache = self.__class__.info_cache
        if self.header_data in info_cache:
            return info_cache[self.header_data][0]

        header_data = self._header_data
        dc = create_info_dataclass(header_data)
        info_cache[header_data] = (dc, create_lazy_class(dc))
        assert issubclass(info_cache[header_data][1], dc)
        return info_cache[header_data][0]
        # return dc

    # @classmethod
    # def adjust_dataclass(cls, dataclass, header):
    #     if not header:
    #         return None
    #     fields = ((field.name, field.type) if field.name != 'info' else ('info', cls._make_info_dataclass(header))
    #               for field in dataclasses.fields(dataclass))
    #     return make_dataclass(fields, dataclass.__name__)


class VCFMatrixBuffer(VCFBuffer):
    dataclass = VCFGenotypeEntry
    # genotype_dataclass = VCFGenotypeEntry
    genotype_encoding = GenotypeRowEncoding

    @classmethod
    def __read_header(cls, file_object):
        # NOT IMPLEMENTED version of read_header()
        prev_line = None
        for line in file_object:
            line = line.decode()
            if line[0] != cls.COMMENT:
                file_object.seek(-len(line), 1)
                break
            prev_line = line
        if prev_line is None:
            return []

        sample_names = prev_line.split("\t")[9:]
        return sample_names

    def __get_data(self):
        data = super().get_data()
        # genotype_data = self._buffer_extractor.get_fixed_length_field(slice(9, None), 3)
        # genotypes = EncodedArray(self.genotype_encoding.encode(genotype_data), self.genotype_encoding)
        # return self.genotype_dataclass(*data.shallow_tuple(), genotypes)

    def _get_field_by_number(self, field_nr: int, field_type: type = object):
        if field_nr != 8:
            assert field_nr < 8, (field_nr, field_type)
            return super()._get_field_by_number(field_nr, field_type)
        genotype_data = self._buffer_extractor.get_fixed_length_field(slice(9, None), 3)
        genotypes = EncodedArray(self.genotype_encoding.encode(genotype_data), self.genotype_encoding)
        return genotypes


# @bnpdataclass
class Genotypes:
    genotype: np.ndarray
    phased: bool


class PhasedVCFMatrixBuffer(VCFMatrixBuffer):
    dataclass = PhasedVCFGenotypeEntry
    lazy_dataclass = create_lazy_class(dataclass)
    genotype_encoding = PhasedGenotypeRowEncoding


class PhasedHaplotypeVCFMatrixBuffer(VCFMatrixBuffer):
    """Encodes genotype info using one column per haplotype (not genotype)"""
    dataclass = PhasedVCFHaplotypeEntry
    lazy_dataclass = create_lazy_class(dataclass)
    genotype_encoding = PhasedHaplotypeRowEncoding
