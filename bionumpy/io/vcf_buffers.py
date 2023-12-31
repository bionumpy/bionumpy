import dataclasses
import logging
from typing import List, Optional

import numpy as np
from npstructures import RaggedArray

from .dump_csv import dump_csv
from .file_buffers import TextBufferExtractor
from .named_text_buffer import NamedBufferExtractor
from .vcf_header import parse_header
from ..bnpdataclass.bnpdataclass import narrow_type
from ..bnpdataclass.lazybnpdataclass import create_lazy_class, ItemGetter, LazyBNPDataClass
from ..encoded_array import EncodedArray, as_encoded_array
from ..bnpdataclass import BNPDataClass, make_dataclass
from ..datatypes import VCFEntry, VCFGenotypeEntry, PhasedVCFGenotypeEntry, PhasedVCFHaplotypeEntry, \
    VCFEntryWithGenotypes, VCFWithInfoAsStringEntry
from ..encodings.vcf_encoding import GenotypeRowEncoding, PhasedGenotypeRowEncoding, PhasedHaplotypeRowEncoding
from .delimited_buffers import DelimitedBuffer
from ..string_array import StringArray

logger = logging.getLogger(__name__)


class GenotypeBufferExtractor(TextBufferExtractor):
    def get_genotypes(self):
        indices = self._starts[..., None] + np.arange(3)
        return self._data[indices]


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


def translate_field_type(info_dict):
    t = info_dict['Type']
    number = info_dict['Number']
    is_list = (number is None) or (number > 1)
    if t == Optional[int] and is_list:
        return List[int]
    elif t == Optional[float] and is_list:
        return List[float]
    elif is_list:
        return str
    return t


def create_info_dataclass(header_data):
    if not header_data:
        return str
    header = parse_header(header_data)
    is_list = lambda val: (val['Number'] is None) or (val['Number'] > 1)
    is_int_list = lambda val: (val['Type'] == Optional[int]) and is_list(val)
    info_fields = [(key, translate_field_type(val)) for key, val in header.INFO.items()]
    dc = make_dataclass(info_fields, "InfoDataclass")
    return dc


class VCFBuffer(DelimitedBuffer):
    '''
    https://samtools.github.io/hts-specs/VCFv4.2.pdf
    '''
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
        elif field_nr == 8:
            return self._extract_genotypes()
        elif field_nr == 9:
            return self._extract_genotype_data()
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
        field_nr = 7
        if (not self._header_data) or ('##INFO' not in self._header_data):
            logger.warning('No header data found or INFO tag missing in header. Cannot parse INFO field.'
                           ' Returning as string. Please use VCFWithInfoAsAstringBuffer to ensure that info field is consistently parsed as string.')
            return self._buffer_extractor.get_field_by_number(field_nr)
        return self._get_dataclass_field(field_nr, self.info_dataclass, self._lazy_info_class)
        # return create_lazy_class(dataclass)(item_getter)

    def _get_dataclass_field(self, field_nr, dataclass, lazy_dataclass):
        text = self._buffer_extractor.get_field_by_number(field_nr, keep_sep=True)
        flat_text = text.ravel()
        delimiters = np.flatnonzero(flat_text == ';') + 1
        offsets = np.insert(np.cumsum(text.lengths), 0, 0)
        all_delimiters = np.sort(np.concatenate([delimiters, offsets]), kind='mergesort')
        delimiter_offsets = np.searchsorted(all_delimiters, offsets)
        dl_lens = np.diff(delimiter_offsets)
        starts = RaggedArray(all_delimiters[:-1].copy(), dl_lens)
        ends = RaggedArray(all_delimiters[1:], dl_lens)
        # starts[:, :-1] = starts[:, :-1] + 1
        ends = ends - 1
        # [:, :-1] = ends[:, :-1] - 1
        assert len(ends)==0 or ends[-1, -1] < len(flat_text), (ends, flat_text)
        lens = ends - starts
        extractor = NamedBufferExtractor(
            flat_text,
            starts,
            lens,
            [f.name for f in dataclasses.fields(dataclass)])
        buf = InfoBuffer(extractor, dataclass)
        item_getter = ItemGetter(buf, dataclass)
        return lazy_dataclass(item_getter)

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

    def _extract_genotypes(self):
        if self._buffer_extractor.n_fields < 9:
            return np.empty((len(self._buffer_extractor), 0), dtype='S')
        byte_array = self._buffer_extractor.get_padded_field(slice(9, None), stop_at=':').raw()
        n_bytes = byte_array.shape[-1]
        if n_bytes == 0:
            bytes = np.empty((len(self._buffer_extractor), self._buffer_extractor.n_fields-9), dtype='S')
        else:
            bytes = byte_array.view(f'>S{n_bytes}').reshape(byte_array.shape[:-1])
        return StringArray(bytes)

    def _extract_genotype_data(self):
        pass

    def get_column_range_as_text(self, col_start, col_end, keep_sep=False):
        if col_start != 8:
            return super().get_column_range_as_text(col_start, col_end, keep_sep=keep_sep)
        return self._buffer_extractor.get_fields_by_range(from_nr=8, to_nr=None, keep_sep=keep_sep)

    @classmethod
    def make_header(cls, data):
        header = ""
        if data.has_context("header"):
            header = data.get_context("header")
        else:
            header='\n'.join([
                '##fileformat=VCFv4.1',
                '\t'.join('#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT'.split())])+ '\n'
        return bytes(header, "ascii")

    @classmethod
    def process_field_for_write(cls, field_name, value):
        if field_name == 'position':
            return value+1
        return super().process_field_for_write(field_name, value)


class VCFBuffer2(VCFBuffer):
    dataclass = VCFEntryWithGenotypes
    lazy_dataclass = create_lazy_class(dataclass)

    @classmethod
    def from_data(cls, data: BNPDataClass) -> "DelimitedBuffer":
        if isinstance(data, LazyBNPDataClass):
            return cls.from_data(data.get_data_object())
        data = dataclasses.replace(data, position=data.position + 1)
        data_dict = [(field.type, getattr(data, field.name)) for field in dataclasses.fields(data)]
        data_dict = data_dict[:-1] + [(str, as_encoded_array(['GT']*len(data)))] + [data_dict[-1]]
        return dump_csv(data_dict, cls.DELIMITER)



class VCFWithInfoAsStringBuffer(VCFBuffer):
    dataclass = VCFWithInfoAsStringEntry


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

class VCFHaplotypeBuffer(VCFBuffer2):
    pass