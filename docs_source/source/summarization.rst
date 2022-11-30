Summarization
~~~~~~~~~~~~~
A key application of BioNumPy is to extract features from sequence datasets. A large set of interesting features can be computed as functions from sequences to scalar values. Examples are kmer-hashing (kmer->hash-value), minimizers(window->hash-value), string/motif-matching (sequence->bool), Position Weight Matrix scores (sequence->float). BioNumPy provides functionality to apply such functions to rolling windows across large sequence sets, through the `RollableFunction` class. By specifying a broadcastable function in the `__call__` method, the `rolling_window` method will apply the function to all windows in a sequence set. Take the `PositionWeightMatrix` class for instance::


    class PositionWeightMatrix(RollableFunction):
        def __init__(self, matrix, encoding=ACTGEncoding):
            self._matrix = np.asanyarray(matrix)
            self.window_size = self._matrix.shape[-1]
            self._indices = np.arange(self.window_size)
            self._encoding = ACTGEncoding

        def __call__(self, sequence: Sequence) -> float:
            sequence = as_sequence_array(sequence, self._encoding)
            scores = self._matrix[sequence, self._indices]
            return scores.sum(axis=-1)

It's `__call__` method specifies how to calculate the score of a sequence. Calling it's rolling_window function will calculate the scores for all windows in a dataset::

    >>> import numpy as np
    >>> sequences = bnp.as_sequence_array(["acgttgta", "gcttca", "gttattc"], encoding=bnp.encodings.ACTGEncoding)
    >>>
    >>> matrix = np.log([[0.1, 0.2],
    ...                  [0.2, 0.3],
    ...                  [0.4, 0.1],
    ...                  [0.3, 0.4]])
    >>> pwm = bnp.position_weight_matrix.PositionWeightMatrix(matrix)
    >>> pwm("ac")
    -3.506557897319982
    >>> pwm(["ac", "cg"])
    array([-3.5065579 , -2.52572864])
    >>> pwm.rolling_window(sequences)
    RaggedArray([[-3.506557897319982, -2.525728644308255, -3.506557897319982, -3.2188758248682006, -1.83258146374831, -3.506557897319982, -2.525728644308255], [-2.4079456086518722, -3.9120230054281455, -3.2188758248682006, -2.120263536200091, -3.2188758248682006], [-3.506557897319982, -3.2188758248682006, -2.525728644308255, -4.605170185988091, -3.2188758248682006, -2.120263536200091]])

Further processing can be done with numpy functions, for instance finding the max score for each sequence in the set::

    >>> pwm.rolling_window(sequences).max(axis=-1)
    array([-1.83258146, -2.12026354, -2.12026354])

