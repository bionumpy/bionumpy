- cite for raggedarray: Awkward Array, Cython, Numba -> Check it out. AA uses cython. Doesn't support some things
- cite/mention biopython for io NCBI retrieval -> How to combine  bionumpiy with biopython
- explain raggedarray implementation (second paragraph). Performance cost -> benchmark, indexing, reductions, elementwise operations, accumulations: see asv
    - asv Ivar lead (1 day maybe pair)
- Pip->pip (typo)

Third paragraph

- Mention supported data file formats/ fasta and genbank files, maybe mention biopython -> What is genbank files (maybe support)

Fourth paragraph

- "In" before "Figure 1" (typo)
- axis labels in figure 1 (typo)

General

- array API: https://data-apis.org/array-api/2022.12/API_specification/index.html (1 day Knut)
    - Two points. Support as much as possible. Support in backend as much as possible (changing backend)
- interplay with torch, jax, numpy -> (1 day pair)
    - change backend to jax (maybe torch), mention cupy.
    - Make work with test suite
    - Hacks in npstructures not in bionumpy
- Autodiff -> (implied)
    - If jax works everything works
- Data types. string/unicode numeric SI?
    -- ACII for bulk data
    -- Encoded arrays for DNA, RNA
    -- Maybe accept unicode

Benchmarks
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
  - Use in npstructures
  - Use also a bit in bionumpy
- y-axis label on each subplot
- Repeats for benchmarks, include variance
    - Make multiple repetitions, and do variance
    - Explain that we benchmark against non-python code
- Make enformer assesment reproducible (if possible)
- "PiP" -> "pip"

Implementation Details

- Figure 2: 2 x 8 = 24 bytes?
- explain vectorization. reduce-at
  - ties in to performance cost earlier
- ufuncs on encoded counts
 -- Maybe implement, low priority
- master -> main

Support alignment
- Data structures [MFA, MAF, CLUSTAL] -> Should be easy
- Algorithms [Wrapper for bwa-mem, minimap2, clustalw, clustalo] -> Only if it's relatively easy (1 day)
    mapper = bwa_mem('hg38.fa')
    aligned_reads = mapper.align('reads.fq')

- Improve documentation -> Timebox (1/2 day together)
    - Add docstrings for all public things
    - Maybe fix typing + add card for it
    - Go a bit through the main documentation and clean up
    - Have 4 concrete stories that should be possible to do with bionumpy
        - interplay with jax and so on
        - end to end pipeline with wrappers
        - immune classification
        - kage
        - (deeptools example)
        - signal plot and matrix



- Blog post -> Chakri (1 day)
-