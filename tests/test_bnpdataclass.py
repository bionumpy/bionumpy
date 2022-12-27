import dataclasses
import pytest
from bionumpy import AminoAcidEncoding, DNAEncoding
from bionumpy.bnpdataclass import bnpdataclass
from bionumpy.bnpdataclass.bnpdataclass import make_dataclass, BNPDataClass
from bionumpy.bnpdataclass.bnpdataclassfunction import bnpdataclassfunction
from numpy.testing import assert_equal
from bionumpy.util.testing import assert_bnpdataclass_equal
import pandas as pd


@bnpdataclass
class Person:
    name: str
    age: int


def test_add_fields():
    @bnpdataclass
    class BaseDC:
        sequence_aa: AminoAcidEncoding

    base_obj = BaseDC(['ACD', "EEA"])
    for field_map_dict in [{"sequence": DNAEncoding}, None]:
        res_obj = base_obj.add_fields({"sequence": ['AA', 'ACT']}, field_map_dict)

        print(res_obj)

        assert all(field.name in ['sequence', 'sequence_aa'] for field in dataclasses.fields(res_obj))
        assert res_obj.sequence.tolist() == ["AA", "ACT"]  # TODO: fix type hinting for fully dynamic stuff


def test_extend():
    @bnpdataclass
    class BaseDC:
        sequence_aa: AminoAcidEncoding

    extended_class = BaseDC.extend((('sequence', DNAEncoding), ('s1', int)))
    assert issubclass(extended_class, BaseDC)
    assert extended_class.__name__ == "DynamicBaseDC"
    assert all(field.name in ['sequence', 'sequence_aa', 's1'] for field in dataclasses.fields(extended_class))


def test_make_dataclass():

    new_cls = make_dataclass([("sequence", DNAEncoding), ('signal1', int)])

    assert issubclass(new_cls, BNPDataClass)
    assert new_cls.__name__ == "DynamicDC"
    assert all(field.name in ['sequence', 'signal1'] for field in dataclasses.fields(new_cls))


def add(a, b):
    return a+b


def test_from_pandas():
    persons = {'name': ['knut', 'marit'],
               'age': [35, 30]}
    df = pd.DataFrame(persons)
    obj = Person(df.name, df.age)
    assert_bnpdataclass_equal(
        obj,
        Person(persons['name'], persons['age']))


@bnpdataclass
class MyClass:
    a: int
    b: int


def test_keyword_init():
    MyClass(a=[10, 20], b=[100, 200])


@pytest.mark.skip("not implemented")
def test_bnpdataclassfunction():
    bnp_add = bnpdataclassfunction("a", "b", (add))
    assert_equal(bnp_add(MyClass([10], [20])), [30])
