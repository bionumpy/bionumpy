import dataclasses
from .bnpdataclass import bnpdataclass, BNPDataClass, make_dataclass

def replace(obj, **kwargs):
    if hasattr(obj, '__replace__'):
        return obj.__replace__(**kwargs)
    return dataclasses.replace(obj, **kwargs)

__all__ = ['bnpdataclass', "BNPDataClass", "make_dataclass"]
