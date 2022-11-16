import scipy.stats

from bionumpy import EncodedRaggedArray, count_encoded
from bionumpy.sequence import get_kmers
from sklearn.decomposition import NMF
from scipy.optimize import minimize
import numpy as np


def estimate_reference_proportions(reference_sequences: EncodedRaggedArray, query_sequences: EncodedRaggedArray):
    reference_kmers = get_kmers(reference_sequences, k=3)
    query_kmers = get_kmers(query_sequences, k=3)
    reference_counts = count_encoded(reference_kmers)
    query_counts = count_encoded(query_kmers, axis=None)
    loss_func = get_loss_func(reference_counts.counts, query_counts.counts)
    print(minimize(loss_func, x0=np.ones(len(reference_sequences))))
    return minimize(loss_func, x0=np.ones(len(reference_sequences))).x ** 2


# def loss(k, rate):
#     (k / rate) - 1


def get_loss_func(reference_counts, query_counts):
    def loss(X):
        return np.sum((X[:, np.newaxis] ** 2 * reference_counts - query_counts) ** 2)

    def poisson_loss(X):
        pmf = scipy.stats.poisson.logpmf(query_counts, np.sum(X[:, np.newaxis] ** 2 * reference_counts, axis=0))
        # print(pmf)
        # print("X ---",X)
        return -np.sum(pmf)

    return poisson_loss
