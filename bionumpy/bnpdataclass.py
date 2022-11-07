import dataclasses
from typing import List
from npstructures.npdataclasses import npdataclass
from npstructures import RaggedArray
import numpy as np
from .encoded_array import EncodedArray, as_encoded_array, EncodedRaggedArray
from .encodings import Encoding, NumericEncoding
from .util import is_subclass_or_instance


def bnpdataclass(base_class: type) -> npdataclass:
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
                assert isinstance(pre_val, (str, list, EncodedArray, EncodedRaggedArray, RaggedArray)), (field, pre_val)
                val = as_encoded_array(pre_val)
            elif is_subclass_or_instance(field.type, Encoding):
                if is_subclass_or_instance(field.type, NumericEncoding):
                    assert isinstance(pre_val, (str, list, EncodedArray, EncodedRaggedArray, RaggedArray, np.ndarray)), (field, pre_val)
                else:
                    assert isinstance(pre_val, (str, list, EncodedArray, EncodedRaggedArray)), (field, pre_val)
                val = as_encoded_array(pre_val, field.type)
            elif field.type == List[int] or field.type == List[bool]:
                if not isinstance(pre_val, RaggedArray):
                    val = RaggedArray(pre_val)
                else:
                    val = pre_val
            else:
                assert False, field.type
            setattr(obj, field.name, val)

    new_class = npdataclass(base_class)
    new_class._implicit_format_conversion = classmethod(_implicit_format_conversion)
    return new_class
