from bionumpy.datatypes import PairsEntry
from .delimited_buffers import DelimitedBuffer


class PairsBuffer(DelimitedBuffer):
    dataclass = PairsEntry
