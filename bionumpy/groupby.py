from npstructures import RaggedView, RaggedArray
import dataclasses
import itertools
import numpy as np
from .npdataclassstream import streamable
from .chromosome_provider import GroupedStream


def get_changes(array):
    if isinstance(array, RaggedArray):
        return get_ragged_changes(array)
    array = array.reshape(len(array), -1)
    return np.flatnonzero(np.all(array[1:]!=array[-1], axis=-1))+1


def get_ragged_changes(ragged_array):
    lengths = ragged_array.shape.lengths
    changes = lengths[1:] != lengths[:-1]
    data = ragged_array.ravel()
    indices =  np.arange(data.size-lengths[-1])
    new_indices = RaggedArray(indices, lengths[:-1])+lengths[:-1, np.newaxis]
    new_indices = np.minimum(new_indices, data.size-1)
    next_row_values = RaggedArray(data[new_indices.ravel()], new_indices.shape)
    eq = next_row_values!=ragged_array[:-1]
    eq = np.any(eq, axis=-1)
    changes |= eq
    return np.flatnonzero(changes)+1


def join_groupbys(grouped_generator):
    double_grouped = itertools.groupby(itertools.chain.from_iterable(grouped_generator), lambda x: x[0])
    return GroupedStream(
        (key, np.concatenate([g[1] for g in groups]))
        for key, groups in double_grouped)


def key_func(x):
    if hasattr(x, "to_string"):
        return x.to_string()
    return str(x)


@streamable(join_groupbys)
def groupby(data, column=None, key=key_func):
    if column is not None:
        assert hasattr(data, column), (data.__class__, dataclasses.fields(data), column)
        keys = getattr(data, column)
    else:
        keys = data
    if np.all(keys[-1] == keys[0]):
        return GroupedStream((key(keys[start]), data[start:])
                             for start in [0])

    changes = get_changes(keys)
    changes = np.append(np.insert(changes, 0, 0), len(data))
    assert np.all(np.diff(changes)>0), changes
    return GroupedStream((key(keys[start]), data[start:end])
                         for start, end in zip(changes[:-1], changes[1:]))
