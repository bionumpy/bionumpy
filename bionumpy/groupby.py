from npstructures import RaggedView, RaggedArray
import itertools
import numpy as np
from .npdataclassstream import streamable
from .chromosome_provider import GroupedStream

def get_changes(array):
    if isinstance(array, RaggedArray):
        return get_ragged_changes(array)
    array = array.reshape(len(array), -1)
    return np.flatnonzero(np.all(array[1:]!=array[-1], axis=-1))

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

    # for key, groups in 
    #     groups = list(groups)
    #     if len(groups) == 1:
    #         yield key, groups[0][1]
    #     else:
    #         yield key, np.concatenate([g[1] for g  in groups])

@streamable(join_groupbys)
def groupby(data, column=None, key=lambda x: x.to_string()):
    if column is not None:
        ragged_array=getattr(data, column)
    else:
        ragged_array=data
    changes = get_changes(ragged_array)
    changes = np.append(np.insert(changes, 0, 0), len(data))
    return GroupedStream((key(ragged_array[start]), data[start:end])
                         for start, end in zip(changes[:-1], changes[1:]))

    # start = 0
    # for change in changes:
    #     yield (str(ragged_array[start]), data[start:change])
    #     start = change
    # yield (str(ragged_array[start]), data[start:len(data)])

