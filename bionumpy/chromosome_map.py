from .chromosome_provider import ChromosomeDictProvider, ChromosomeStreamProvider
import logging
logger = logging.getLogger(__name__)



class ChromosomeMap:
    def __init__(self, reduction=None):
        self._reduction = reduction

    def __call__(self, func):
        def mapped(*args, **kwargs):
            stream_indices = [i for i, a in enumerate(args) if isinstance(a, ChromosomeStreamProvider)]
            dict_indices = [i for i, a in enumerate(args) if isinstance(a, ChromosomeDictProvider)]
            stream_keys = [key for key, val in kwargs.items() if isinstance(val, ChromosomeStreamProvider)]
            dict_keys = [key for key, val in kwargs.items() if isinstance(val, ChromosomeDictProvider)]
            assert len(stream_indices)+len(stream_keys) <= 1
    
            if len(stream_indices) == 1:
                stream = args[stream_indices[0]]
            elif len(stream_keys) == 1:
                stream = kwargs[stream_keys[0]]
            new_args = list(args)
            dicts = [args[i] for i in dict_indices]
            dicts_kw = [kwargs[key] for key in dict_keys]
            for chromosome, data in stream:
                logger.info(f"Running '{func.__name__}' on chromosome {chromosome}")
                for i, d in zip(dict_indices, dicts):
                    new_args[i] = d[chromosome]
    
                for key, d in zip(dict_keys, dicts_kw):
                    kwargs[key] = d[chromosome]
                for i in stream_indices:
                    new_args[i] = data
                for key in stream_keys:
                    kwargs[key] = data
                    
                yield func(*new_args, **kwargs)
                
        if self._reduction is None:
            return mapped
        return lambda *args, **kwargs: self._reduction(mapped(*args, **kwargs))


def sorted_streams_juggler(streams):
    """TODO"""
    n_streams = len(streams)
    cur_values = [next(stream, None) for stream in streams]
    while True:
        next_thing = min(cur_values)
