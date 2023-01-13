from typing import Union

from .stream import BnpStream, NpDataclassStream
from . import groupby
from ..bnpdataclass import BNPDataClass
import logging
import sys

from ..datatypes import ChromosomeSize

logger = logging.getLogger(__name__)


def human_key_func(chrom_name):
    assert chrom_name.startswith("chr"), chrom_name
    parts = chrom_name[3:].split("_", maxsplit=1)
    assert len(parts) <= 2, chrom_name
    is_numeric = 1-parts[0].isdigit()
    b = parts[0] if is_numeric else int(parts[0])
    c = parts[-1] if len(parts) == 2 else ""
    return (is_numeric, b, c)


def sort_dict_by_key(dictionary, key=None):
    return {key: dictionary[key] for key in sorted(dictionary.keys())}


class StreamError(Exception):
    pass


class SequenceSizes(dict):
    pass


class SynchedStream(BnpStream):
    """
    Groups a stream based on a grouping attribute (defaults to "chromosome"), 
    and returns data in the order specified in `contig_order`, with empty datastructures
    if the stream has no attribute for that contig. The column of the grouping attribute
    needs to be sorted with the same order as `contig_order`. 
    
    If it's not, and you're only working with one stream, it's recomended to change the 
    order of `contig_order` rather than sorting the data file. The sorting mehtod is
    usually either normal string sort, or sorted using `alpha_numeric_key_func`.
    """
    def __init__(self, stream, contig_order):
        self._stream = stream
        self._contig_order = contig_order
        self._grouping_attribute = "chromosome"
        self._has_default = True
        self._default_value = stream.dataclass.empty()
        self._key_func = lambda x: x

    def set_grouping_attribute(self, attribute_name):
        self._grouping_attribute = attribute_name

    def set_key_function(self, key_function):
        self._key_func = key_function

    def __iter__(self):
        grouped = groupby(self._stream, self._grouping_attribute)
        cur_contig_idx = 0
        used_names = []
        seen_contig_names = set([])
        for name, data in grouped:
            logger.debug(f"handling data for {name}")
            sys.stdout.flush()
            sys.stderr.flush()
            name = self._key_func(name)
            if name in seen_contig_names:
                raise StreamError(f"Sort order discrepancy between stream and contig. {name} already occured in {seen_contig_names}")
            used_names.append(name)
            if name not in self._contig_order:
                raise StreamError(f"Stream had value not present in contig order: {name} ({self._contig_order})")
            while cur_contig_idx < len(self._contig_order) and (name != self._contig_order[cur_contig_idx]):
                if self._has_default:
                    logger.info(f"Data for contig {self._contig_order[cur_contig_idx]} missing; using default value: {self._default_value}")
                    sys.stdout.flush()
                    sys.stderr.flush()
                    yield self._default_value
                    seen_contig_names.add(self._contig_order[cur_contig_idx])
                    cur_contig_idx += 1
                else:
                    raise StreamError(
                        f"Stream's next element ({name}) is not {self._contig_order[cur_contig_idx]} as expected by contig order. And no default value is set for the stream")
            if name == self._contig_order[cur_contig_idx]:
                yield data
                seen_contig_names.add(self._contig_order[cur_contig_idx])
                cur_contig_idx += 1
            else:
                raise StreamError(f"Stream element {name} not present in contig order : {self._contig_order}")
        if cur_contig_idx < len(self._contig_order):
            if self._has_default:

                for i in range(cur_contig_idx, len(self._contig_order)):
                    logger.info(f"Data for contig {self._contig_order[i]} missing")
                    sys.stdout.flush()
                    sys.stderr.flush()
                    yield self._default_value
            else:
                raise StreamError(
                    f"Stream's next element ({name}) is not {self._contig_order[cur_contig_idx]} as expected by contig order. And no default value is set for the stream")
        remains = next(grouped, None)
        if remains is not None:
            raise StreamError(f"Stream element {remains[0]} not present in contig ordere {self._contig_order}")
        
    def __repr__(self):
        with_default = f" with default_value {self._default_value}" if self._has_default else ""
        return f"SynchedStream over {self._contig_order}{with_default}.\n{self._stream}"
            
    __str__ = __repr__
    def set_default(self, default_value):
        self._has_default = True
        self._default_value = default_value


class IndexedStream(BnpStream):
    """Creates a stream from a dict-like object that generates
    values from the dict in the order of `contig_order`
    """
    def __init__(self, lookup, contig_order):
        self._lookup = lookup
        self._contig_order = contig_order

    def __iter__(self):
        return (self._lookup[name] for name in self._contig_order)

    def __repr__(self):
        return f"IndexedStream over contigs: {self._contig_order}\n{self._lookup}"

    __str__ = __repr__


class MultiStream:
    """ Class to handle multiple streams/data sources that works on the same
    set of reference sequences

    Examples
    --------
    >>> from bionumpy.streams import NpDataclassStream
    >>> from bionumpy.datatypes import Interval
    >>> positions = {"chr1": [1, 2, 3], "chr2": [3, 4, 5]}
    >>> sequence_sizes = {"chr1": 10, "chr2": 15}
    >>> stream = NpDataclassStream([Interval(["chr1"]*2, [0, 4], [5, 10]), Interval(["chr2"], [2], [13])], dataclass=Interval)
    >>> multistream = MultiStream(sequence_sizes, positions=positions, intervals=stream)
    >>> for pos, interval in zip(multistream.positions, multistream.intervals):
    ...     print(pos)
    ...     print(interval)
    [1, 2, 3]
    Interval with 2 entries
                   chromosome                    start                     stop
                         chr1                        0                        5
                         chr1                        4                       10
    [3, 4, 5]
    Interval with 1 entries
                   chromosome                    start                     stop
                         chr2                        2                       13
    """
    def __init__(self, sequence_sizes: Union[SequenceSizes, ChromosomeSize], **kwargs):
        """Synch different data streams with data on the same reference
    
        Take streams or dict-like objects and return a stream that gives
        you all the data for one contig/chromosome at the time
    
        Parameters
        ----------
        sequence_sizes : SequenceSizes
            A mapping from contig-name to contig-size
        **kwargs : key, value pairs for each data source
    
        """
        if isinstance(sequence_sizes, dict):
            sequence_names = list(sequence_sizes.keys())
            sequence_lengths = list(sequence_sizes.values())
        elif isinstance(sequence_sizes, ChromosomeSize):
            sequence_names = sequence_sizes.name.tolist()
            sequence_lengths = sequence_sizes.size.tolist()
        else:
            raise Exception("MultiStream needs sequence_sizes as either a "
                            "dict-like object or a ChromosomeSize object, not %s" % sequence_sizes)

        self._streams = {}
        self.lengths = BnpStream(sequence_lengths)
        self.sequence_names = BnpStream(sequence_names)
        for keyword, value in kwargs.items():
            if isinstance(value, BNPDataClass):
                value = NpDataclassStream([value], value.__class__)
            if isinstance(value, BnpStream):
                self.__dict__[keyword] = SynchedStream(value, sequence_names)
            elif hasattr(value, "__getitem__"):
                self.__dict__[keyword] = IndexedStream(value, sequence_names)
            else:
                raise ValueError(f"Only streams and dict-like objects can be used in MultiStream. {keyword}:{value}")

    def set_defaults(self, **kwargs):
        for keyword, default_value in kwargs.items():
            assert keyword in self.__dict__
            self.__dict__[keyword].set_default(default_value)

    def set_key_functions(self, **kwargs):
        for keyword, key_function in kwargs.items():
            assert keyword in self.__dict__
            self.__dict__[keyword].set_key_function(key_function)

    @staticmethod
    def human_key_func(chrom_name):
        assert chrom_name.startswith("chr"), chrom_name
        parts = chrom_name[3:].split("_", maxsplit=1)
        assert len(parts) <= 2, chrom_name
        is_numeric = 1-parts[0].isdigit()
        b = parts[0] if is_numeric else int(parts[0])
        c = parts[-1] if len(parts) == 2 else ""
        return (is_numeric, b, c)

    @staticmethod
    def sort_dict_by_key(dictionary, key=None):
        return {key: dictionary[key] for key in sorted(dictionary.keys(), key=key)}
