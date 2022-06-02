import numpy as np
a = np.array([[1, 2, 3],
              [4, 5, 6]])
b = np.array([[1, 10, 100],
              [2, 20, 200]])
a+b

c = np.array([2, 3, 5])
a+c

d = np.array([10, 20])
a+d

d = np.array([[10], [20]])
a+d


v = np.array([2, 3, 5])
w = np.array([3, 5, 7])
v.dot(w)
(v*w).sum()

v.dot([1])

big = np.arange(30).reshape(2, 5, 3)
v.dot(big)


def count_gs(sequence):
    return np.sum(sequence == "G")


def one_hot(letter, alphabet):
    return letter == alphabet

def one_hot(letter, alphabet):
    return letter[:, np.newaxis] == alphabet
