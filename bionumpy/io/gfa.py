from typing import List
from npstructures import RaggedArray
from ..bnpdataclass import bnpdataclass
from .delimited_buffers import DelimitedBuffer
from .strops import split, str_to_int
import numpy as np


@bnpdataclass
class GfaPath:
    name: str
    node_ids: List[int]
    directions: List[int]


class GfaPathBuffer(DelimitedBuffer):

    def get_data(self):
        name = self.get_text(1, fixed_length=False)
        nodes_lists = self.get_text(2, keep_sep=True, fixed_length=False)
        nodes_lists[:, -1] = ","
        lengths = np.sum(nodes_lists == ",", axis=-1)
        all_node_texts = split(nodes_lists.ravel()[:-1], ",")
        int_text = all_node_texts[:, :-1]
        node_ids = str_to_int(int_text)
        directions = np.where(all_node_texts[:, -1]=="+", 1, -1)
        node_ids = RaggedArray(node_ids, lengths)
        directions = RaggedArray(directions, lengths)
        data =  GfaPath(name, node_ids, directions)
        return data
