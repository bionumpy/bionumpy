import logging
from .stream import BnpStream

logger = logging.getLogger(__name__)


def grouped_dict(attribute_name=None):
    def decorator(base_class):
        base_class.grouped_dict_attribute = attribute_name
        return base_class
    return decorator


class grouped_stream(BnpStream):
    def __init__(self, stream, attribute_name=None):
        self.attribute_name = attribute_name
        super().__init__(stream)


class chromosome_map:
    """
    Decorator that applies a the function to the data of each chromosome
    if a ChromosomeProvider instance is given as an argument.

    If no ChromosomeProvider instance is given as argument, the behaviour of the
    function is unchanged
    """

    def __init__(self, reduction=None):
        self._reduction = reduction

    def get_args(
        self, args, kwargs, stream_indices, dict_indices, stream_keys, dict_keys
    ):
        if len(stream_indices) == 1:
            stream = args[stream_indices[0]]
            self.attribute_name = stream.attribute_name
        elif len(stream_keys) == 1:
            stream = kwargs[stream_keys[0]]
            self.attribute_name = stream.attribute_name
        elif len(dict_indices) > 0:
            stream_indices.append(dict_indices.pop(0))
            self.attribute_name = args[stream_indices[-1]].grouped_dict_attribute
            stream = args[stream_indices[-1]].items()
        elif len(dict_keys) > 0:
            stream_keys.append(dict_keys.pop(0))
            self.attribute_name = args[stream_keys[-1]].grouped_dict_attribute
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

    @staticmethod
    def is_grouped_dict(obj):
        return hasattr(obj, "grouped_dict_attribute")

    def __call__(self, func):
        def log(sequence):
            for chromosome, args, kwargs in sequence:
                logger.info(f"Running {func.__name__} on chromosome {chromosome}")
                yield chromosome, args, kwargs

        def mapped(*args, **kwargs):
            stream_indices = [
                i for i, a in enumerate(args) if isinstance(a, grouped_stream)
            ]
            dict_indices = [
                i for i, a in enumerate(args) if self.is_grouped_dict(a)
            ]
            stream_keys = [
                key
                for key, val in kwargs.items()
                if isinstance(val, grouped_stream)
            ]
            dict_keys = [
                key
                for key, val in kwargs.items()
                if self.is_grouped_dict(val)
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
                ret = grouped_stream(((chromosome, func(*args, **kwargs))
                                      for chromosome, args, kwargs in log(new_args)))
            elif is_dict:
                ret = grouped_dict()(dict((chromosome, func(*args, **kwargs))
                                          for chromosome, args, kwargs in log(new_args)))
            else:
                assert False
            if self._reduction is None:
                return ret
            return self._reduction(
                (m[1] for m in ret))
        return mapped
