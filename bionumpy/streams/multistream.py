from .stream import BnpStream
from ..groupby import groupby


class StreamError(Exception):
    pass


class SequenceSizes(dict):
    pass


class SynchedStream(BnpStream):
    def __init__(self, stream, contig_order):
        self._stream = stream
        self._contig_order = contig_order
        self._grouping_attribute = "chromosome"
        self._has_default = False

    def set_grouping_attribute(self, attribute_name):
        self._grouping_attribute = attribute_name

    def __iter__(self):
        grouped = groupby(self._stream, self._grouping_attribute)
        cur_contig_idx = 0
        for name, data in grouped:
            if name not in self._contig_order:
                raise StreamError(f"Stream had value not present in contig order: {name} ({self._contig_order})")
            if name == self._contig_order[cur_contig_idx]:
                yield data
                cur_contig_idx += 1
            else:
                if self._has_default:
                    yield self._default_value
                else:
                    raise StreamError(
                        f"Stream's next element is not {self._contig_order[cur_contig_idx]} as expected by contig order. And no default value is set for the stream")
        if cur_contig_idx < len(self._contig_order):
            if self._has_default:
                for _ in range(cur_contig_idx, len(self._contig_order)):
                    yield self._default_value
            else:
                raise StreamError(
                    f"Stream's next element is not {self._contig_order[cur_contig_idx]} as expected by contig order. And no default value is set for the stream")

    def __repr__(self):
        with_default = f" with default_value {self._default_value}" if self._has_default else ""
        return f"SynchedStream over {self._contig_order}{with_default}.\n{self._stream}"
            
    def set_default(self, default_value):
        self._has_default = True
        self._default_value = default_value


class IndexedStream(BnpStream):
    def __init__(self, lookup, contig_order):
        self._lookup = lookup
        self._contig_order = contig_order

    def __iter__(self):
        return (self._lookup[name] for name in self._contig_order)

    def __repr__(self):
        return f"IndexedStream over contigs: {self._contig_order}\n{self._lookup}"


class MultiStream:
    def __init__(self, sequence_sizes: SequenceSizes, **kwargs):
        """Synch different data streams with data on the same reference
    
        Take streams or dict-like objects and return a stream that gives
        you all the data for one contig/chromosome at the time
    
        Parameters
        ----------
        sequence_sizes : SequenceSizes
            A mapping from contig-name to contig-size
        **kwargs : key, value pairs for each data source
    
        Examples
        --------
        FIXME: Add docs.
        """
        self._streams = {}
        self.lengths = list(sequence_sizes.values())
        for keyword, value in kwargs.items():
            if isinstance(value, BnpStream):
                self.__dict__[keyword] = SynchedStream(value, list(sequence_sizes.keys()))
            elif hasattr(value, "__getitem__"):
                self.__dict__[keyword] = IndexedStream(value, list(sequence_sizes.keys()))
            else:
                raise ValueError(f"Only streams and dict-like objects can be used in MultiStream. {keyword}:{value}")

    def set_defaults(self, **kwargs):
        for keyword, default_value in kwargs.items():
            assert keyword in self.__dict__
            self.__dict__[keyword].set_default(default_value)
