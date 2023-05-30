import npstructures as nps
from .. import EncodedRaggedArray, EncodedArray


def ragged_slice(array: EncodedRaggedArray, starts=None, ends=None) -> EncodedRaggedArray:
    slice = nps.ragged_slice(array.ravel(), starts, ends)
    return EncodedRaggedArray(EncodedArray(slice.ravel(), array.encoding), slice.shape)