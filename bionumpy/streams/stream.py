class BnpStream:
    def __init__(self, stream, first_buffer=None):
        self._stream = iter(stream)
        self._first_buffer = first_buffer
        self._opened = False

    @property
    def is_opened(self):
        return self._opened

    def __iter__(self):
        return self

    def __next__(self):
        self._opened = True
        return next(self._stream)

    def _peek(self):
        if self._peeked:
            return self._buffer
        self._buffer = next(self._stream, None)
        self._peeked = True
        return self._buffer

    def __str__(self):
        status = "opened" if self._opened else "unopened"
        return f"""\
{status.capitalize()} bionumpy-stream of data-buffers. Next buffer:
{self._peek()}"""


class ArrayStream(BnpStream):
    def __init__(self, stream, first_buffer=None):
        self._stream = stream
        self._first_buffer = first_buffer

    def ravel(self):
        return ArrayStream(arr.ravel() for arr in self)


class NpDataclassStream(BnpStream):
    """
    Class to hold a stream/generator of npdataclass objects.
    Works with streamable/bnp_broadcast so that functions decorated
    with those decorator will apply the function to each element in
    a NpDataclassStream. To concatenate all the chunks in the generator
    into one npdataclass, use `NpDataclassStream.join()`
    """
    def __init__(self, stream, dataclass):
        self._stream = stream
        self._opened = False
        self._peeked = False
        self._buffer = None
        self._dataclass = dataclass

    def _peek(self):
        if self._peeked:
            return self._buffer
        self._buffer = next(self._stream, None)
        self._peeked = True
        return self._buffer


        return next(self._stream)

    def __iter__(self):
        return self
    
    def __next__(self):
        return next(self._stream)

        self._opened = True
        if self._peeked:
            self._peeked = False
            buf = self._buffer
            self._buffer = None
            return buf

    def __str__(self):
        status = "opened" if self._opened else "unopened"
        return f"""\
{status.capitalize()} stream of {self._dataclass} buffers. Next buffer:
{self._peek()}"""

    def __getattr__(self, attribute_name):
        if attribute_name not in {f.name for f in dataclasses.fields(self._dataclass)}:
            raise Exception(f"{self._dataclass} has no attribute {attribute_name}")
        return ArrayStream(getattr(chunk, attribute_name) for chunk in self)
