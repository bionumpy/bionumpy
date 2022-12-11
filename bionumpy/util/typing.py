from typing import List, Union, NewType
from ..encoded_array import EncodedArray, EncodedRaggedArray

SingleEncodedArrayLike = NewType("SingleEncodedArrayLike", Union[str, EncodedArray])
EncodedRaggedArrayLike = NewType("EncodedRaggedArrayLike", Union[List[str], EncodedRaggedArray, List[EncodedArray]])
EncodedArrayLike = NewType("EncodedArrayLike", Union[SingleEncodedArrayLike, EncodedRaggedArrayLike])
