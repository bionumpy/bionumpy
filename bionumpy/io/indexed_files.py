from pathlib import PurePath
import os
from .indexed_fasta import IndexedFasta, create_index
from .files import bnp_open
from .delimited_buffers import DelimitedBuffer
from .multiline_buffer import FastaIdx


class IndexBuffer(DelimitedBuffer):
    sep = "\t"
    dataclass = FastaIdx


def open_indexed(filename: str) -> IndexedFasta:
    """Open an indexed fasta (for now) file with random access

    If an index is not already present for the file, create it

    Parameters
    ----------
    filename : str
        The filename of the file

    Returns
    -------
    IndexedFasta
        An Indexed fasta object that supports random access on
        chromosome or intervals

    Examples
    --------
    >>> from bionumpy import open_indexed
    >>> reference = open_indexed("example_data/small_genome.fa")
    >>> reference
    Indexed Fasta File with chromosome sizes: {'0': 80, '1': 80, '2': 80, '3': 80}
    >>> reference["1"]
    encoded_array('gcttggtatgaaaacccatc...')
    >>> from bionumpy.datatypes import Interval
    >>> intervals = Interval.from_entry_tuples([("1", 10, 20), ("2", 20, 30)])
    >>> reference.get_interval_sequences(intervals)
    encoded_ragged_array(['aaaacccatc',
                          'ggccgttttt'])
    """

    path = PurePath(filename)
    suffix = path.suffixes[-1]
    index_file_name = path.with_suffix(path.suffix + ".fai")
    assert suffix in (".fa", ".fasta"), "Only fasta supported for indexed read"
    if not os.path.isfile(index_file_name):
        bnp_open(index_file_name, "w", buffer_type=IndexBuffer).write(create_index(path))
    return IndexedFasta(filename)
