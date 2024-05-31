import npstructures as nps
from .. import EncodedRaggedArray, EncodedArray


def ragged_slice(array: EncodedRaggedArray, starts=None, ends=None) -> EncodedRaggedArray:
    """
    Slice a ragged array column-wise.
    Parameters
    ----------
    array : EncodedRaggedArray
        The array to slice

    starts : np.ndarray
        The start indices of the slices
    ends : np.ndarray, optional
        The end indices of the slices. If not provided, the slices will be taken to the end of the array.

    Returns
    -------
    EncodedRaggedArray
        The sliced array

    Examples
    --------
    >>> import numpy as np
    >>> import bionumpy as bnp
    >>> seqs = bnp.as_encoded_array(["ACGT", "ACGTT"])
    >>> starts = np.array([0, 1])
    >>> ends = np.array([2, 3])
    >>> sliced = bnp.ragged_slice(seqs, starts, ends)
    >>> print(sliced)
    AC
    CG
    """
    sliced_data = nps.ragged_slice(array.ravel(), starts, ends)
    return EncodedRaggedArray(EncodedArray(sliced_data.ravel(), array.encoding), sliced_data.shape)