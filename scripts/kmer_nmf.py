import scipy.stats
import bionumpy as bnp
from bionumpy.encoded_array import EncodedRaggedArray, count_encoded, as_encoded_array
from bionumpy.dna import reverse_compliment
from bionumpy.sequence import get_kmers
from bionumpy.simulate import rnaseq
from sklearn.decomposition import NMF
from scipy.optimize import minimize
import numpy as np
import plotly.express as px
from bionumpy.cli import run_as_commandline


def estimate_reference_proportions(reference_sequences: EncodedRaggedArray, query_sequences: EncodedRaggedArray):
    reference_kmers = get_kmers(reference_sequences, k=3)
    query_kmers = get_kmers(query_sequences, k=3)
    reference_counts = count_encoded(reference_kmers) + count_encoded(get_kmers(reverse_compliment(reference_sequences), k=3))
    query_counts = count_encoded(query_kmers, axis=None)
    print("reference_counts", reference_counts.counts)
    print("query_counts", query_counts.counts)
    loss_func = get_loss_func(reference_counts.counts, query_counts.counts)
    return minimize(loss_func, x0=np.ones(len(reference_sequences))).x ** 2 * 2 * (reference_sequences.shape.lengths - 2)


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


def test():
    transcriptome = as_encoded_array(["ACGT", "GGGATC", "AAATCG", "TTTTACG"], bnp.DNAEncoding)
    rnaseq_reads = rnaseq.simulate_rnaseq(transcriptome,
                                          rnaseq.RNASeqSimulationSettings(fragment_size=4, read_length=3,
                                                                          transcript_counts=[9, 18, 27, 36]))
    estimated_reads = estimate_reference_proportions(transcriptome, rnaseq_reads)
    return estimated_reads

xx|x|z|
def main(re):
    

if __name__ == '__main__':
    run_as_commandline(estimate_reference_proportions(
    est_reads = test()
    fig = px.scatter(x=est_reads, y=[9, 18, 27, 36])
    fig.show()
