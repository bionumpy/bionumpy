from .chromosome_provider import ChromosomeDictProvider, ChromosomeStreamProvider

def chromosome_map(reduction=None):
    def decorator(func):
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
                print(f"Running chromosome {chromosome}")
                for i, d in zip(dict_indices, dicts):
                    new_args[i] = d[chromosome]
    
                for key, d in zip(dict_keys, dicts_kw):
                    kwargs[key] = d[chromosome]
                for i in stream_indices:
                    new_args[i] = data
                for key in stream_keys:
                    kwargs[key] = data
                    
                yield func(*new_args, **kwargs)
                
        if reduction is None:
            return mapped
        return lambda *args, **kwargs: reduction(mapped(*args, **kwargs))
    return decorator


def sorted_streams_juggler(streams):
    """TODO"""
    n_streams = len(streams)
    cur_values = [next(stream, None) for stream in streams]
    while True:
        next_thing = min(cur_values)
