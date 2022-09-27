Rollable Functions
==================

In many cases, we want to take a simple function mapping a sequence to some scalar, and map it to all subsequences of a given length in a set of sequences. Examples of this is:
* Hashing kmers
* Use a position weight matrix to compute a score for a sequence
* Find all occurances of a substring in a sequence set

bioNumpy provides this functionality throught the `RollableFunction` class. All you have to to is subclass the `RollableFunction` class, and write a broadcastable version of the sequence function as the `__call__` method. A call to the `rolling_window` method will then apply the function to all the subsequences of length `window_size` in the sequence set. `window_size` can either be set as `self.window_size` or passed as argument to the `rolling_window` method.

For instance, if we want to check for instances of "CGGT" in a set of sequences, we can use the following::

    from bionumpy.rollable import RollableFunction
    from bionumpy.sequences import as_sequence_array
    import numpy as np
    
    
    class StringMatcher(RollableFunction):
        def __init__(self, matching_sequence):
            self._matching_sequence = as_sequence_array(matching_sequence)
    
        def __call__(self, sequence):
            return np.all(sequence == self._matching_sequence, axis=-1)

The `__call__` function here just checks that all the letters in the sequence are equal to the corresponding letters in the matching sequence. Specifying `axis=-1` for the all function makes the function broadcastable::

    >>> matcher = StringMatcher("CGGT")
    >>> matcher("CGGT")
    Sequence(True)

Giving a sequence of different length to the `__call__` function returns `False`, since the sequneces are then not equal::

    >>> matcher("CGGTA")
    <stdin>:7: DeprecationWarning: elementwise comparison failed; this will raise an error in the future.
    False

However we can use the `rolling_window` method to match every subsequence of length 4 to "CGGT"::

    >>> matcher.rolling_window("CGGTA")
    array([ True, False])
    >>> matcher.rolling_window(["CGGTA", "ACGGTG"])
    RaggedArray([[True, False], [False, True, False]])

For examples of rollable function implementations see:
* `Minimizers`
* `KmerEncoding`
* `PositionWeightMatrix`
