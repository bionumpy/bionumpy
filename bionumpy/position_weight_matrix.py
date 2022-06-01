from .rollable import RollableFunction

LobProbabilities = None

class PositionWeightMatrix:
    def __init__(self, matrix, encoding=ACGTEncoding):
        self._matrix = np.asanyarray(matrix)
        self.window_size = self._matrix.shape[-1]
        self._indices = np.arange(self.window_size)
        self._alphabet = alphabet

    def __call__(self, sequence):
        scores = self._matrix[sequence, self._indices]
        return scores.sum(axis=-1)

    @property
    def range(self):
        return LogProbabilites

    def sample_domain(self, n):
        return self.encoding.sample_domain(n*self.window_size).reshape(-1, self.window_size)



@rolling_window_function
def rolling_pwm(sequence, window_size, pwm):
    assert window_size == pwm.window_size
    return pwm.calculate_score(sequence)
