import numpy as np


def parse_jaspar_line(line):
    letter, rest = line.split(maxsplit=1)
    rest = rest.strip()[1:-1].split()
    counts = [int(n) for n in rest]
    return letter.strip(), counts


def read_jaspar_matrix(filename):
    f = open(filename)
    _ = f.readline()
    alphabet, matrix = zip(*(parse_jaspar_line(line) for line in f))
    alphabet = "".join(alphabet)
    matrix = np.array(matrix, dtype="float")
    return alphabet, matrix


if __name__ == "__main__":
    filename = "/home/knut/Downloads/MA0001.1.jaspar"
    alphabet, matrix = read_jaspar_matrix(filename)
    pwm = get_pwm
