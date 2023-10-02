import numpy as np


class IntegerEncoding:
    def __init__(self, nan_symbol: str = None):
        self._nan_symbol = nan_symbol

    def _encode(self, data: np.ndarray):
        encoded = data-ord('0')
        if self._nan_symbol is not None:
            encoded[data==ord(self._nan_symbol)] = np.nan
        return encoded

    def _decode(self, data: np.ndarray):
        return
