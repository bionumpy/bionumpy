from npstructures import RaggedView, RaggedArray
import numpy as np

def get_changes(ragged_array):
    print(ragged_array)
    lengths = ragged_array.shape.lengths
    changes = lengths[1:] != lengths[:-1]
    data = ragged_array.ravel()
    indices =  np.arange(data.size-lengths[-1])
    new_indices = RaggedArray(indices, lengths[:-1])+lengths[:-1, np.newaxis]
    new_indices = np.minimum(new_indices, data.size-1)
    print(new_indices)
    next_row_values = RaggedArray(data[new_indices.ravel()], new_indices.shape)
    print(next_row_values)
    eq = next_row_values!=ragged_array[:-1]
    eq = np.any(eq, axis=-1)
    changes |= eq
    return np.flatnonzero(changes)+1
    
