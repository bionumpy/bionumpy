from .chromosome_provider import (
    ChromosomeDictProvider,
    ChromosomeStreamProvider,
    PureChromosomeDictProvider,
)
import logging

logger = logging.getLogger(__name__)


class ChromosomeMap:
    def __init__(self, reduction=None):
        self._reduction = reduction

    def get_args(
        self, args, kwargs, stream_indices, dict_indices, stream_keys, dict_keys
    ):
        if len(stream_indices) == 1:
            stream = args[stream_indices[0]]
        elif len(stream_keys) == 1:
            stream = kwargs[stream_keys[0]]
        elif len(dict_indices) > 0:
            stream_indices.append(dict_indices.pop(0))
            stream = args[stream_indices[-1]].items()
        elif len(dict_keys) > 0:
            stream_keys.append(dict_keys.pop(0))
            stream = kwargs[stream_keys[-1]].items()

        new_args = list(args)
        dicts = [args[i] for i in dict_indices]
        dicts_kw = [kwargs[key] for key in dict_keys]
        for chromosome, data in stream:
            for i, d in zip(dict_indices, dicts):
                new_args[i] = d[chromosome]
            for key, d in zip(dict_keys, dicts_kw):
                kwargs[key] = d[chromosome]
            for i in stream_indices:
                new_args[i] = data
            for key in stream_keys:
                kwargs[key] = data

            yield chromosome, new_args, kwargs

    def __call__(self, func):
        def log(sequence):
            for chromosome, args, kwargs in sequence:
                logger.info(f"Running {func.__name__} on chromosome {chromosome}")
                yield chromosome, args, kwargs

        def mapped(*args, **kwargs):
            stream_indices = [
                i for i, a in enumerate(args) if isinstance(a, ChromosomeStreamProvider)
            ]
            dict_indices = [
                i for i, a in enumerate(args) if isinstance(a, ChromosomeDictProvider)
            ]
            stream_keys = [
                key
                for key, val in kwargs.items()
                if isinstance(val, ChromosomeStreamProvider)
            ]
            dict_keys = [
                key
                for key, val in kwargs.items()
                if isinstance(val, ChromosomeDictProvider)
            ]
            is_stream = len(stream_indices) + len(stream_keys) > 0
            is_dict = not is_stream and (len(dict_keys) + len(dict_indices) > 0)

            if not (is_stream or is_dict):
                return func(*args, **kwargs)
            assert len(stream_indices) + len(stream_keys) <= 1
            new_args = self.get_args(
                args, kwargs, stream_indices, dict_indices, stream_keys, dict_keys
            )
            if is_stream:
                return ChromosomeStreamProvider(
                    (chromosome, func(*args, **kwargs))
                    for chromosome, args, kwargs in log(new_args)
                )
            # self.get_args(args, kwargs, stream_indices, dict_indices, stream_keys, dict_keys)))
            elif is_dict:
                return PureChromosomeDictProvider(
                    (chromosome, func(*args, **kwargs))
                    for chromosome, args, kwargs in log(new_args)
                )
            # self.get_args(args, kwargs, stream_indices, dict_indices, stream_keys, dict_keys))
            assert False

        if self._reduction is None:
            return mapped
        return lambda *args, **kwargs: self._reduction(
            (m[1] for m in mapped(*args, **kwargs))
        )


def sorted_streams_juggler(streams):
    """TODO"""
    n_streams = len(streams)
    cur_values = [next(stream, None) for stream in streams]
    while True:
        next_thing = min(cur_values)
