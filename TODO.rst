- cite for raggedarray: Awkward Array, Cython, Numba
- cite/mention biopython for io NCBI retrieval
- explain raggedarray implementation (second paragraph). Performance cost
- Pip->pip

Third paragraph

- Mention supported data file formats/ fasta and genbank files, maybe mention biopython

Fouth paragraph

- "In" before "Figure 1"
- axis labels in figure 1

General

- array API: https://data-apis.org/array-api/2022.12/API_specification/index.html
  - from importlib.metadata import entry_points
  - return hasattr(x, '__array_namespace__')
  - https://data-apis.org/array-api/2022.12/verification_test_suite.html
  - doesn't support list inputs
  - mutation is not allowed
  - access device in raggedarrays
  - dlpack
  - .devce, device=None, to_device
  - in-place operations?
  - we can't do T, mT
  - __array_namespace__
  - __index__
  - __int__
  - broadcasting - we only support 2 dims
  - array creation algorithms
  - positional only?
  - must return zero-dim array on indexing (this differs from numpy)
  - boolean indexing is voluntary
  - take
  - manipulation functions
  - argmin, argmax, nonzero, where
  - set functions
  - argsort, sort
  - max, mean, min, prod, std, sum, var, all, any
  - check type promotion, should be delegated entirely to numpy
  - __array_api_version__
  - extension name
  
- interplay with torch jax, numpy
- Autodiff
- Data types. string/unicode numeric SI?- Supporting Info
- Benchmarks
- Use asv: https://asv.readthedocs.io/en/v0.6.1/)
- y-axis label on each subplot
- Repeats for benchmarks, include variance
- Make enformer assesment reproducible
- "PiP" -> "pip"

Implementation Details

- Figure 2: 2 x 8 = 24 bytes?
- explain vectorization. reduce-at
- ufuncs on encoded counts
- masteer -> main                                      

Support alignment
- Data structures
- Algorithms ??

- Impprove documentation
