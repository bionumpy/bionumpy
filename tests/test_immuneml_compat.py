import dataclasses
from typing import Any, Dict

import pytest
import bionumpy as bnp
from bionumpy import AminoAcidEncoding, DNAEncoding
from bionumpy.bnpdataclass import bnpdataclass

text='''\
sequence_aa\tsequence\tv_call
MMM\tACGT\tV1
CCC\tACGT\tV2
'''


@pytest.fixture
def header(type_dict):
    return '\t'.join(type_dict.keys()) + '\n'


@pytest.fixture()
def full_text(header):
    text = '''\
a\tb\tc\t1\t2\ta\ta\ta\tacgt\taaaa\tid\tv\n\
'''
    return header + text


@pytest.fixture
def file_name(tmp_path):
    name = tmp_path / 'tmp1234.tsv'
    with open(name, 'w') as f:
        f.write(text)
    return name

@pytest.fixture
def empty_file_name(header, tmp_path):
    name = tmp_path / 'empty.csv'
    with open(name, 'w') as f:
        f.write(header)
    return name

@pytest.fixture
def full_file_name(full_text, tmp_path):
    name = tmp_path / 'tmp1234full.tsv'
    with open(name, 'w') as f:
        f.write(full_text)
    return name

@bnpdataclass
class ReceptorSequence:
    sequence_aa: AminoAcidEncoding
    sequence: DNAEncoding
    v_call: str


def raise_func(*args, **kwargs):
    return True


@pytest.fixture
def type_dict():
    return {'cell_id': str, 'chain': str, 'cmv': str,
            'coeliac': bool, 'duplicate_count': int,
            'frame_type': str, 'j_call': str,
            'region_type': str,
            'sequence': DNAEncoding, 'sequence_aa': AminoAcidEncoding,
            'sequence_id': str, 'v_call': str}

@pytest.fixture
def simple_type_dict():
    return {'sequence': DNAEncoding, 'v_call': str, 'sequence_aa': AminoAcidEncoding}


def make_dynamic_seq_set_dataclass(type_dict: Dict[str, Any]):
    dc = dataclasses.make_dataclass('DynamicSequenceSet', fields=type_dict.items(),
                                    namespace={
                                        'get_row_by_index': raise_func,
                                        'get_single_row_value': raise_func,
                                        'to_dict': raise_func,
                                        'get_rows_by_indices': raise_func})
    return bnpdataclass(dc)


def create_buffer_type_from_field_dict(type_dict: Dict[str, Any]) -> bnp.io.delimited_buffers.DelimitedBuffer:
        dataclass = make_dynamic_seq_set_dataclass(type_dict)
        return bnp.io.delimited_buffers.get_bufferclass_for_datatype(
            dataclass, delimiter='\t', has_header=True)





@pytest.fixture()
def buffer_type(simple_type_dict):
    return create_buffer_type_from_field_dict(simple_type_dict)


@pytest.fixture()
def full_buffer_type(type_dict):
    return create_buffer_type_from_field_dict(type_dict)


def test_read(file_name):
    buffer_type = bnp.io.delimited_buffers.get_bufferclass_for_datatype(
        ReceptorSequence, delimiter='\t',
        has_header=True)
    with bnp.open(str(file_name), buffer_type=buffer_type) as file:
        obj = file.read()
    print(obj.get_data_object())
    print(obj.sequence_aa,
          obj.sequence, obj.v_call)


def test_read2(file_name, buffer_type):
    with bnp.open(str(file_name), buffer_type=buffer_type) as file:
        obj = file.read()
    print(obj.get_data_object())
    print(obj.sequence_aa, obj.sequence, obj.v_call)


def test_read2(full_buffer_type, full_file_name):
    with bnp.open(str(full_file_name), buffer_type=full_buffer_type) as file:
        obj = file.read()
    print(obj.tolist())
    print(obj.get_data_object())
    print(obj.sequence_aa, obj.sequence, obj.v_call)

def test_has_methods(full_buffer_type, full_file_name):
    import bionumpy.config
    bionumpy.config.LAZY = False
    with bnp.open(str(full_file_name), buffer_type=full_buffer_type) as file:
        obj = file.read()

    print(obj.get_row_by_index(0))
    bionumpy.config.LAZY = True

def test_read_empty(empty_file_name, buffer_type):
    data = bnp.open(str(empty_file_name), buffer_type=buffer_type).read()
    assert len(data) == 0
