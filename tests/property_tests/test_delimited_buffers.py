from bionumpy.bnpdataclass import bnpdataclass
from npstructures.testing import assert_npdataclass_equal
from bionumpy.io.delimited_buffers import get_bufferclass_for_datatype
from typing import List
from .strategies import integers, floats, ascii_text
from hypothesis import strategies as st
from hypothesis import given, example
from functools import partial
import dataclasses


type_to_strategy = {int: integers,
                    str: partial(ascii_text, min_size=1),
                    float: floats,
                    List[int]: lambda: st.lists(integers, min_size=1)}


@bnpdataclass
class MyDataclass:
    name: str
    age: int


def table_strategies(dataclass):
    fixed_dict = {field.name: type_to_strategy[field.type]() for field in dataclasses.fields(dataclass)}
    return st.lists(st.fixed_dictionaries(fixed_dict), min_size=1)

def table_to_dataclass(dataclass, table):
    return dataclass(*[
        [row[field.name] for row in table]
        for field in dataclasses.fields(dataclass)])


@given(table_strategies(MyDataclass))
@example(tables=[{'age': 0, 'name': '0'}, {'age': -1, 'name': '0'}])
def test_to_from_data(tables):
    data = table_to_dataclass(MyDataclass, tables)
    buffer_class = get_bufferclass_for_datatype(MyDataclass)
    buf = buffer_class.from_data(data)
    new_data = buffer_class.from_raw_buffer(buf).get_data()
    assert_npdataclass_equal(new_data, data)
