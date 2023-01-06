import numpy as np
from ..sequence.position_weight_matrix import PWM


def parse_jaspar_line(line):
    letter, rest = line.split(maxsplit=1)
    rest = rest.strip()[1:-1].split()
    counts = [int(n) for n in rest]
    return letter.strip(), counts


def read_jaspar_matrix(filename):
    f = open(filename)
    _ = f.readline()
    pwm = dict((parse_jaspar_line(line) for line in f))
    return PWM.from_dict(pwm)
    # alphabet, matrix = zip(*(parse_jaspar_line(line) for line in f))
    # 
    # alphabet = "".join(alphabet)
    # matrix = np.array(matrix, dtype="float")
    # return PWM.from_dict({char: row for
    # return alphabet, matrix
