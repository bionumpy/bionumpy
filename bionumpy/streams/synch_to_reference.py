from .stream import BnpStream
from ..groupby import groupby


class SequenceSizes:
    pass


def synch_to_reference(sequence_sizes: SequenceSizes, **kwargs):
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
    chromosome_names = sequence_sizes.keys()
    return_values = {}
    for key, value in sequence_sizes.items():
        if isinstance(value, BnpStream):
            new_value = groupby(value, "chromosome")
        elif hasattr(value, "__getitem__"):
            pass
