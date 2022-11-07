class BnpStream:
    def __init__(self, stream):
        self._stream = iter(stream)
        self._next_buffer = next(self._stream, None)
        self._opened = False

    @property
    def is_opened(self):
        return self._opened

    def __iter__(self):
        return self

    def __next__(self):
        self._opened = True
        if self._next_buffer is None:
            raise StopIteration
        result = self._next_buffer
        self._next_buffer = next(self._stream, None)
        return result

    def __str__(self):
        status = "opened" if self._opened else "unopened"
        return f"""\
{status.capitalize()} bionumpy-stream of data-buffers. Next buffer:
{self._next_buffer}"""

    def __repr__(self):
        return f"self.__class__.__name__({self._next_buffer}, {Ellipsis})"


class ChunkStream(BnpStream):
    pass


class ArrayStream(ChunkStream):
    pass


class NpDataclassStream(BnpStream):
    """
    Class to hold a stream/generator of npdataclass objects.
    Works with streamable/bnp_broadcast so that functions decorated
    with those decorator will apply the function to each element in
    a NpDataclassStream. To concatenate all the chunks in the generator
    into one npdataclass, use `NpDataclassStream.join()`
    """
    def __init__(self, stream, dataclass=None):
        super().__init__(stream)
        self.dataclass = dataclass

    def __getattr__(self, attribute_name):
        return ArrayStream(getattr(chunk, attribute_name) for chunk in self)
