LAZY = True
STRING_ARRAY = True


class ConfigContext:
    def __init__(self, **kwargs):
        self._kwargs = kwargs

    def __enter__(self):
        global LAZY, STRING_ARRAY
        self._old_kwargs = {k: globals()[k] for k in self._kwargs}
        globals().update(self._kwargs)

    def __exit__(self, exc_type, exc_value, traceback):
        global LAZY, STRING_ARRAY
        globals().update(self._old_kwargs)
