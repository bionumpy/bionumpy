import functools

from . import bnp_open
from .bnpdataclass.bnpdataclass import BNPDataClass
import inspect
try:
    import typer
except ImportError:
    typer = None


def check_arguments(function):
    arguments = inspect.getfullargspec(function)
    print(arguments)
    new_annotations = {name: str if isinstance(t, BNPDataClass) else t for name, t in arguments.annotations.items()}
    if BNPDataClass(arguments.annotations['return']):
        pass


class CliWrapper:
    '''Convert all arguments that are typed with BNPDataclass into filename options'''

    def __init__(self, *args, **kwargs):
        self._args = args
        self._kwargs = kwargs

    def __call__(self, function):
        argspec = inspect.getfullargspec(function)
        do_write = 'return' in argspec.annotations and issubclass(argspec.annotations['return'], BNPDataClass)

        def is_bnpdataclass(name: int) -> bool:
            return issubclass(argspec.annotations[name], BNPDataClass)

        @functools.wraps(function)
        def new_func(*args, **kwargs):
            new_args = [bnp_open(arg).read() if is_bnpdataclass(argspec.args[i]) else arg for i, arg in enumerate(args)]
            new_kwargs = {k: bnp_open(v).read() if is_bnpdataclass(k) else v for k, v in kwargs.items() if k != 'output'}
            return_val = function(*new_args, **new_kwargs)
            if do_write:
                bnp_open(kwargs['output'], "w").write(return_val)
        sig = inspect.signature(function)
        new_parameters = [val.replace(annotation=str) if issubclass(val.annotation, BNPDataClass) else val for key, val in sig.parameters.items()]
        if do_write:
            new_parameters.append(inspect.Parameter('output', inspect.Parameter.KEYWORD_ONLY, annotation=str, default=None))
        new_sig = sig.replace(parameters=new_parameters, return_annotation=sig.empty)
        new_func.__signature__ = new_sig

        # Add the new annotations to the function
        annotations = {name: str if is_bnpdataclass(name) else t for name, t in argspec.annotations.items() if name != 'return'}
        if do_write:
            annotations['output'] = str
        new_func.__annotations__ = annotations
        return new_func
