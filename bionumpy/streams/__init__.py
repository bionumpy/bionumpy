from .stream import BnpStream
from .decorators import streamable
NpDataclassStream = None

mean, bincount, histogram, quantile = (None for _ in range(4))

__all__ = ["BnpStream", "streamable", "decorators",
           "mean", "bincount", "histogram", "quantile"]
