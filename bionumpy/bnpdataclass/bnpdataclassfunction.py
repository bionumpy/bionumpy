from .bnpdataclass import BNPDataClass
import dataclasses


class bnpdataclassfunction:
    def __init__(self, *args):
        arg_names = args

    def __call__(self, func):
        def new_func(data_object, *args, **kwargs):
            pass


def apply_to_npdataclass(attribute_name):
    def decorator(func):
        def new_func(np_dataclass, *args, **kwargs):
            if not isinstance(np_dataclass, BNPDataClass):
                return func(np_dataclass, *args, **kwargs)
            result = func(getattr(np_dataclass, attribute_name), *args, **kwargs)
            return dataclasses.replace(np_dataclass, **{attribute_name: result})
        return new_func
    return decorator


