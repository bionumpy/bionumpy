from typing import List
from .encoded_array import EncodedArray, EncodedRaggedArray

SingleEncodedArrayLike = str | EncodedArray
EncodedRaggedArrayLike = List[str] | EncodedRaggedArray | List[EncodedArray]
EncodedArrayLike = SingleEncodedArrayLike | EncodedRaggedArrayLike
