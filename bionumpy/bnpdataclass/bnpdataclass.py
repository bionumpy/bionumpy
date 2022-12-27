import dataclasses
from typing import List, Type, Dict
from numpy.typing import ArrayLike
from npstructures.npdataclasses import npdataclass, NpDataClass
from npstructures import RaggedArray
import numpy as np
from ..encoded_array import EncodedArray, EncodedRaggedArray
from ..encoded_array import as_encoded_array
from ..encodings import Encoding, NumericEncoding
from ..util import is_subclass_or_instance


class BNPDataClass(NpDataClass):

    @classmethod
    def extend(cls, fields: tuple, name: str = None) -> Type['BNPDataClass']:
        """
        Parameters
        ----------
        fields: tuple
            A tuple in format (field_name, field_type) for the new fields to be added
        name: str
            The optional user-defined name for the new class

        Returns
        --------
        BNPDataClass with added fields

        Examples
        ---------

        >>> from bionumpy.bnpdataclass import bnpdataclass
        >>> from bionumpy.encodings import AminoAcidEncoding, DNAEncoding
        >>> @bnpdataclass
        ... class BaseDC:
        ...     sequence_aa: AminoAcidEncoding

        >>> extended_class = BaseDC.extend((('sequence', DNAEncoding), ('s1', int)))
        >>> assert all(field.name in ['sequence', 'sequence_aa', 's1'] for field in dataclasses.fields(extended_class))
        >>> print([field.name for field in dataclasses.fields(extended_class)])
        ['sequence_aa', 'sequence', 's1']

        """
        cls_name = name if name is not None else f"Dynamic{cls.__name__}" if cls.__name__[:7] != 'Dynamic' else cls.__name__
        return bnpdataclass(dataclasses.make_dataclass(cls_name, bases=(cls,), fields=fields))

    def add_fields(self, fields: Dict[str, ArrayLike], field_type_map: dict = None) -> 'BNPDataClass':
        """
        Parameters
        ----------
        fields: dict
            a dictionary in containing the names of the new fields as keys and lists of values for each of the field as values

        field_type_map: dict
            a dictionary with field names as keys and types as values; for basic types, they can be inferred from the data and don't need to be
            specified; but for fields that need to use some of the encodings, the specific encoding can be provided here

        Returns
        --------
        BNPDataClass object with added fields with the provided values

        Examples
        ---------

        >>> from bionumpy.bnpdataclass import bnpdataclass
        >>> from bionumpy.encodings import AminoAcidEncoding, DNAEncoding
        >>> @bnpdataclass
        ... class BaseDC:
        ...     sequence_aa: AminoAcidEncoding

        >>> base_obj = BaseDC(['ACD', "EEA"])
        >>> res_obj = base_obj.add_fields({"sequence": ['AA', 'ACT']}, field_type_map={'sequence': DNAEncoding})
        >>> print(res_obj)
        DynamicBaseDC with 2 entries
                      sequence_aa                 sequence
                              ACD                       AA
                              EEA                      ACT

        """
        for name in fields.keys():
            if not name.isidentifier():
                raise TypeError(f"Field name must be a valid identifier (No whitespace or dots and such): {name}")
        fields_with_types = _extract_field_types(fields, field_type_map)
        new_class = self.__class__.extend(tuple(fields_with_types.items()))
        return new_class(**{**vars(self), **fields})

    @classmethod
    def from_entry_tuples(cls, tuples):
        return cls(*(list(c) for c in zip(*tuples)))


def bnpdataclass(base_class: type) -> Type[BNPDataClass]:
    """Create a `bnpdataclass` from a class with fields specified

    A wrapper around `@npdataclass` that includes implicit format conversion
    for strings and other encoded data. `@npdataclass` is again a wrapper
    around `dataclasses.dataclass` but where all the fields are assumed to be
    objects that supports advanced numpy indexing so that any dataclass objects
    are also indexible like a numpy array.

    `bnpdataclass` classes are meant to be dataclasses and so should not have 
    any methods other than those implicitly given by `npdataclass`.

    Parameters
    ----------
    base_class : type
        Base class that defines the fields of the dataclass.

    Returns
    -------
    npdataclass
        `bnpdataclass` object that supports numpy like indexing

    Examples
    --------
    >>> from bionumpy.bnpdataclass import bnpdataclass
    >>> @bnpdataclass
    ... class Person:
    ...       name: str
    ...       age: int
    ... 
    >>> data = Person(["Knut", "Ivar", "Geir"], [35, 30, 40])
    >>> print(data)
    Person with 3 entries
                         name                      age
                         Knut                       35
                         Ivar                       30
                         Geir                       40

    >>> print(data[[0,2]])
    Person with 2 entries
                         name                      age
                         Knut                       35
                         Geir                       40
    
    >>> print(data[[False,True, False]])
    Person with 1 entries
                         name                      age
                         Ivar                       30
    """

    class NewClass(npdataclass(base_class), BNPDataClass):
        @classmethod
        def _implicit_format_conversion(cls, obj: npdataclass):
            """Convert the data in given in the init into numpy like data

            This is convenience functionionality that converts e.g. list of strings
            into EncodedRaggedArray objects, and numeric data into `np.ndarray` objects.

            Called by the `__init__` method from `npdataclass`.

            Parameters
            ----------
            cls : 3
                The class this is called from
            obj : npdataclass
                The partially initialzed object from `npdataclass` `__init__`

            """
            for field in dataclasses.fields(obj):
                pre_val = getattr(obj, field.name)
                if field.type in (int, float):
                    val = np.asanyarray(pre_val)
                elif field.type == str:
                    assert isinstance(pre_val, (str, list, EncodedArray, EncodedRaggedArray, RaggedArray, np.ndarray)) or hasattr(pre_val, 'to_numpy'), (field, pre_val, type(pre_val))
                    val = as_encoded_array(pre_val)
                elif is_subclass_or_instance(field.type, Encoding):
                    if is_subclass_or_instance(field.type, NumericEncoding):
                        assert isinstance(pre_val, (str, list, EncodedArray, EncodedRaggedArray, RaggedArray, np.ndarray)), \
                            (field, pre_val, type(pre_val))
                    else:
                        assert isinstance(pre_val, (str, list, EncodedArray, EncodedRaggedArray)), (field, pre_val)
                    # must do as_encoded and not explicit encode as pre_val might already
                    # be encoded
                    val = as_encoded_array(pre_val, field.type)
                elif field.type == List[int] or field.type == List[bool]:
                    if not isinstance(pre_val, RaggedArray):
                        val = RaggedArray(pre_val)
                    else:
                        val = pre_val
                else:
                    assert False, field.type
                setattr(obj, field.name, val)

    NewClass.__name__ = base_class.__name__
    NewClass.__qualname__ = base_class.__qualname__

    return NewClass


def make_dataclass(fields: list, name: str = "DynamicDC") -> Type[BNPDataClass]:
    """
    Constructs a dynamic dataclass from a list of attributes

    Parameters
    ----------
    fields: list
        a list of tuples in format (field_name, field_type) to be used to construct the dynamic bnp dataclass

    name: str
        optional name of new class

    Returns
    -------
    new BNPDataClass

    """
    return bnpdataclass(dataclasses.make_dataclass(name, fields=fields))


def _extract_field_types(fields_with_values: dict, field_type_map: dict = None) -> dict:
    fields = {}
    for field_name in fields_with_values.keys():
        _assert_all_same_type(fields_with_values[field_name])

        if field_type_map is not None and field_name in field_type_map:
            field_type = field_type_map[field_name]
        elif isinstance(fields_with_values[field_name][0], EncodedArray):
            field_type = type(fields_with_values[field_name][0].encoding)
        else:
            field_type = type(fields_with_values[field_name][0])

        if fields_with_values[field_name] is not None:
            fields[field_name] = field_type

    return fields


def _assert_all_same_type(values):
    original_type = type(values[0])
    assert all(isinstance(val, original_type) for val in values), (original_type, [type(val) for val in values])
