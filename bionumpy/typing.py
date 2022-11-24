from typing import List, Union
from .encoded_array import EncodedArray, EncodedRaggedArray

SingleEncodedArrayLike = Union(str, EncodedArray)
EncodedRaggedArrayLike = Union(List[str], EncodedRaggedArray, List[EncodedArray])
EncodedArrayLike = Union(SingleEncodedArrayLike, EncodedRaggedArrayLike)
