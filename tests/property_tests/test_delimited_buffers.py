from bionumpy.bnpdataclass import bnpdataclass
from bionumpy.util.testing import assert_bnpdataclass_equal
# from npstructures.testing import assert_bnpdataclass_equal
from bionumpy.io.delimited_buffers import get_bufferclass_for_datatype
import bionumpy.datatypes as dt
from typing import List
from .strategies import integers, floats, ascii_text
from hypothesis import strategies as st
from hypothesis import given, example
from functools import partial
import dataclasses
import pytest

type_to_strategy = {int: integers,
                    str: partial(ascii_text, min_size=1),
                    float: lambda: floats().filter(lambda x: abs(x) > 10**(-15)),
                    List[int]: partial(st.lists, elements=integers(), min_size=1),
                    List[bool]: partial(st.lists, elements=st.booleans(), min_size=1),
                    }


@bnpdataclass
class MyDataclass:
    name: str
    age: int
    money: float
    child_ages: List[int]
    child_gender: List[bool]


def table_strategies(dataclass):
    fixed_dict = {field.name: type_to_strategy[field.type]() for field in dataclasses.fields(dataclass)}
    return st.lists(st.fixed_dictionaries(fixed_dict), min_size=1)


def table_to_dataclass(dataclass, table):
    return dataclass(*[
        [row[field.name] for row in table]
        for field in dataclasses.fields(dataclass)])



# @given(table_strategies(MyDataclass))
# @example(tables=[{'age': 0, 'name': '0'}, {'age': -1, 'name': '0'}])
# @example(tables=[{'age': 0, 'child_ages': [0], 'money': 0.0, 'name': '0'}])
@pytest.mark.skip("Skipped because requires encoding already encoded array with different encoding. Not supported")
def _test_to_from_data(tables):
    data = table_to_dataclass(MyDataclass, tables)
    buffer_class = get_bufferclass_for_datatype(MyDataclass)
    buf = buffer_class.from_data(data)
    file_buffer = buffer_class.from_raw_buffer(buf)
    print(file_buffer._data)
    new_data = file_buffer.get_data()
    assert_bnpdataclass_equal(new_data, data)


test_to_from_data = given(table_strategies(MyDataclass))(_test_to_from_data)

# for datatype in datatypes:
#     setattr(__main__

# objs = (getattr(dt, name) for name in dir(dt) if not name.startswith("_"))
# datatypes = [obj for obj in objs if hasattr(obj, "shallow_tuple")]

