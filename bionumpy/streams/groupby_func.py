from npstructures import RaggedView, RaggedArray
import dataclasses
import itertools
import numpy as np
from . import streamable, grouped_stream
from ..bnpdataclass import bnpdataclass
from ..encoded_array import EncodedArray
from ..encodings.string_encodings import StringEncoding
from ..string_array import StringArray


def get_changes(array):
    if isinstance(array, EncodedArray) and isinstance(array.encoding, StringEncoding):
        return np.flatnonzero(array.raw()[1:] != array.raw()[:-1])+1
    if isinstance(array, RaggedArray):
        return get_ragged_changes(array)
    elif isinstance(array, StringArray):
        return np.flatnonzero(array.raw()[1:] != array.raw()[:-1])+1
    array = array.reshape(len(array), -1)
    return np.flatnonzero(np.all(array[1:]!=array[:-1], axis=-1))+1


def get_ragged_changes(ragged_array):
    lengths = ragged_array.lengths
    changes = lengths[1:] != lengths[:-1]
    data = ragged_array.ravel()
    indices =  np.arange(data.size-lengths[-1])
    new_indices = RaggedArray(indices, lengths[:-1])+lengths[:-1, np.newaxis]
    new_indices = np.minimum(new_indices, data.size-1)
    next_row_values = RaggedArray(data[new_indices.ravel()], new_indices.shape)
    eq = next_row_values != ragged_array[:-1]
    eq = np.any(eq, axis=-1)
    changes |= eq
    return np.flatnonzero(changes)+1


def join_groupbys(grouped_generator):
    double_grouped = itertools.groupby(itertools.chain.from_iterable(grouped_generator), lambda x: x[0])
    def f(groups):
        groups_ = [g[1] for g in groups]
        return np.concatenate(groups_)

    grouped_ = ((key, f(groups)) for key, groups in double_grouped)
    return grouped_stream(grouped_,
                          grouped_generator.attribute_name if hasattr(grouped_generator, "attribute_name") else None)


def key_func(x):
    if hasattr(x, "to_string"):
        return x.to_string()
    return str(x)


@streamable(join_groupbys)
def groupby(data: bnpdataclass, column: str=None, key: callable = key_func):
    """Group the data according to the values in `column`

    This functions behaves similarily to `itertools.groupby`, but
    requires the values to be sorted beforehand. It will return a
    generator yielding (name, data) pairs, where name is the value of
    the groupby'ed column, and data is a bnpdataclass with all the
    entries for that value. This functions works on both `npdataclass`
    and `NpDataClassStream`.

    The main intended use for this function is to group data per chromosome/contig.

    Parameters
    ----------
    data : bnpdataclass
        The data to be grouped
    column : str
        The name of the attribute that should be used for the grouping
    key : callable
        A function to be called on the grouped by attribute for the returned tuple
    Examples
    --------
    >>> from bionumpy.datatypes import Interval
    >>> intervals = Interval(["1", "1", "2", "2", "3"], [10, 14, 5, 17, 3], [15, 20, 10, 20, 10])
    >>> print(intervals)
    Interval with 5 entries
                   chromosome                    start                      end
                            1                       10                       15
                            1                       14                       20
                            2                        5                       10
                            2                       17                       20
                            3                        3                       10
    >>> for name, data in groupby(intervals, "chromosome"):
    ...       print(name)
    ...       print(data)
    ...
    1
    Interval with 2 entries
                   chromosome                    start                      end
                            1                       10                       15
                            1                       14                       20
    2
    Interval with 2 entries
                   chromosome                    start                      end
                            2                        5                       10
                            2                       17                       20
    3
    Interval with 1 entries
                   chromosome                    start                      end
                            3                        3                       10

    """
    if column is not None:
        assert hasattr(data, column), (data.__class__, data, column)
        keys = getattr(data, column)
    else:
        keys = data
    if (isinstance(keys, EncodedArray) or (hasattr(keys, "lengths") and keys.lengths[-1] == keys.lengths[0])) and np.all(keys[-1] == keys[0]):
        return grouped_stream(((key(keys[start]), data[start:]) for start in [0]), column)
                                       

    changes = get_changes(keys)
    changes = np.append(np.insert(changes, 0, 0), len(data))
    assert np.all(np.diff(changes)>0), changes
    return grouped_stream(((key(keys[start]), data[start:end])
                           for start, end in zip(changes[:-1], changes[1:])),
                          column)
