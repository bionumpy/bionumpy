import dataclasses
from itertools import accumulate
import numpy as np

class NpDataClass:
    def shallow_tuple(self):
        return tuple(getattr(self, field.name) for field in dataclasses.fields(self))

    def __getitem__(self, idx):
        return self.__class__(*[f[idx] for f in self.shallow_tuple()])

    def __len__(self):
        return len(self.shallow_tuple()[0])

    def __array_function__(self, func, types, args, kwargs):
        if func==np.concatenate:
            objects = args[0]
            tuples = [o.shallow_tuple() for o in objects]
            new_tuple = tuple
            return self.__class__(*(np.concatenate(list(t)) for t in zip(*tuples)))
        if func == np.equal:
            one, other = args
            return all(np.equal(s, o) for s, o in zip(one.shallow_tuple(), other.shallow_tuple()))
            
        return NotImplemented

    def __iter__(self):
        return (self.__class__(*comb) for comb in zip(*self.shallow_tuple()))


class VarLenArray:
    def __init__(self, array):
        self.array = array
        self.shape = self.array.shape
        self.dtype = self.array.dtype

    def __array_function__(self, func, types, args, kwargs):
        if func==np.concatenate:
            arrays = [v.array for v in args[0]]
            lens = [len(arg) for arg in arrays]
            sizes = [arg.shape[-1] for arg in arrays]
            max_size = max(sizes)
            if all(size == max_size for size in sizes):
                return self.__class__(np.concatenate(arrays))
            ret = np.zeros((sum(lens), max_size), dtype=self.array.dtype)
            for end, l,  a, size in zip(accumulate(lens), lens, arrays, sizes):
                ret[end-l:end, -size:] = a
            return self.__class__(ret)
        if func == np.equal:
            raise Exception()
        return NotImplemented


    def __eq__(self, other):
        return self.array==other.array

    def __repr__(self):
        return f"VarLen({repr(self.array)})"

    def __neq__(self, other):
        return self.array != other.array

    def __len__(self):
        return len(self.array)

    def __array__(self, *args, **kwargs):
        return self.array
    
    def __iter__(self):
        return iter(self.array)

    def __getitem__(self, idx):
        return self.__class__(self.array[idx])



