import numpy as np

n_letters = 128

def _get_power_array(n, mod):
    l = [1]
    for _ in range(n-1):
        l.append((l[-1]*n_letters % mod))
    return np.array(l)

def get_ascii_hash(encoded_array, mod):
    powers = _get_power_array(encoded_array.shape[-1], mod)
    return np.sum((powers*encoded_array.raw()) % mod) % mod
    
