import dataclasses

import numpy as np

from ..encoded_array import EncodedArray
from ..bnpdataclass import BNPDataClass
from ..datatypes import VCFEntry, VCFGenotypeEntry, PhasedVCFGenotypeEntry, PhasedVCFHaplotypeEntry
from ..encodings.vcf_encoding import GenotypeRowEncoding, PhasedGenotypeRowEncoding, PhasedHaplotypeRowEncoding
from .delimited_buffers import DelimitedBuffer


class VCFBuffer(DelimitedBuffer):
    dataclass = VCFEntry

    def get_data(self):
        data = super().get_data()
        data.position -= 1
        return data

    def get_field_by_number(self, field_nr: int, field_type: type=object):
        val = super().get_field_by_number(field_nr, field_type)
        if field_nr == 1:
            val -= 1
        return val

    @classmethod
    def from_data(cls, data: BNPDataClass) -> "DelimitedBuffer":
        data = dataclasses.replace(data, position=data.position+1)
        return super().from_data(data)


class VCFMatrixBuffer(VCFBuffer):
    dataclass = VCFEntry
    genotype_dataclass = VCFGenotypeEntry
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

    def get_data(self):
        data = super().get_data()
        genotypes = self.get_column_range_as_text(9, self._n_cols, keep_sep=True)
        genotypes = EncodedArray(self.genotype_encoding.encode(genotypes), self.genotype_encoding)
        return self.genotype_dataclass(*data.shallow_tuple(), genotypes)

    def get_field_by_number(self, field_nr: int, field_type: type=object):
        if field_nr != 9:
            return super().get_field_by_number(field_nr, field_type)
        else:
            genotypes = self.get_column_range_as_text(9, self._n_cols, keep_sep=True)
            genotypes = EncodedArray(self.genotype_encoding.encode(genotypes), self.genotype_encoding)
            return genotypes


#@bnpdataclass
class Genotypes:
    genotype: np.ndarray
    phased: bool




class PhasedVCFMatrixBuffer(VCFMatrixBuffer):
    dataclass = VCFEntry
    genotype_dataclass = PhasedVCFGenotypeEntry
    genotype_encoding = PhasedGenotypeRowEncoding


class PhasedHaplotypeVCFMatrixBuffer(VCFMatrixBuffer):
    """Encodes genotype info using one column per haplotype (not genotype)"""
    dataclass = VCFEntry
    genotype_dataclass = PhasedVCFHaplotypeEntry
    genotype_encoding = PhasedHaplotypeRowEncoding
