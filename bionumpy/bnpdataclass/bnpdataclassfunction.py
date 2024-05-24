from .bnpdataclass import BNPDataClass
import dataclasses


class bnpdataclassfunction:
    def __init__(self, *args):
        arg_names = args

    def __call__(self, func):
        def new_func(data_object, *args, **kwargs):
            pass


def replace(obj: BNPDataClass, **kwargs) -> BNPDataClass:
    '''Replace the values of a dataclass with new values

    Parameters
    ----------
    obj : BNPDataClass
        The dataclass to be replaced
    kwargs : dict
        The new values to be replaced

    Returns
    -------
    BNPDataClass
        The new dataclass with the replaced values

    Examples
    --------
    >>> import bionumpy as bnp
    >>> entry = bnp.SequenceEntry(['seq1'], ['acgt'])
    >>> entry
    SequenceEntry with 1 entries
                         name                 sequence
                         seq1                     acgt
    >>> bnp.replace(entry, name=['seq2'])
    SequenceEntry with 1 entries
                         name                 sequence
                         seq2                     acgt
    '''
    if hasattr(obj, '__replace__'):
        return obj.__replace__(**kwargs)
    return dataclasses.replace(obj, **kwargs)


def apply_to_npdataclass(attribute_name):
    def decorator(func):
        def new_func(np_dataclass, *args, **kwargs):
            if not isinstance(np_dataclass, BNPDataClass):
                return func(np_dataclass, *args, **kwargs)
            result = func(getattr(np_dataclass, attribute_name), *args, **kwargs)
            return replace(np_dataclass, **{attribute_name: result})

        return new_func

    return decorator
