from .stream import BnpStream
from typing  import Generator
import numpy as np

def _chunk_entries(stream: BnpStream, n_entries: int) -> Generator:

    b = []
    buffer_size = 0
    for chunk in stream:
        b.append(chunk)
        buffer_size += len(chunk)
        if buffer_size >= n_entries:
            total = np.concatenate(b)
            yield total[:n_entries]
            b = [total[n_entries:]]
            buffer_size = len(b[0])
    if buffer_size:
        yield np.concatenate(b)


def chunk_entries(stream: BnpStream, n_entries: int) -> BnpStream:
    """Chunk a stream into fixed number of entries

    Parameters
    ----------
    stream : BnpStream
    n_entries : int

    Returns
    -------
    BnpStream
    """
    return stream.__class__(_chunk_entries(stream, n_entries))
