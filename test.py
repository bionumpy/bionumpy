
import numpy as np
import cupy as cp

x = cp.arange(10)
print(x)

y = x[cp.asanyarray([0, 2])]
print(y)
