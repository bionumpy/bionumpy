import os
import gzip
import numpy as np
import yaml
from importlib import reload
import bionumpy.sequence.kmers
# from npstructures import RaggedArray
import npstructures.raggedarray as npsa
import npstructures.raggedshape as npss
import bionumpy.encodings
import milem.categorical as EM
# from milem import SparseEM, MultiCategorical
reload(EM)
reload(bionumpy.encodings)
reload(bionumpy.sequence.kmers)
reload(npsa)
reload(npss)
RaggedArray = npsa.RaggedArray


def read_repertoire(data_filename, meta_filename):
    y = yaml.load(open(meta_filename))
    print(y)
    i = y["field_list"].index("sequence_aas")
    d = np.load(gzip.open(data_filename))
    seqs = [row[i] for row in d]
    return RaggedArray([[ord(c) for c in seq] for seq in seqs])


def get_positive_repertoires_filenames(folder):
    return get_labelled_repertoirs_filenames(folder)


def get_labelled_repertoirs_filenames(folder, label=True):
    files = [f for f in os.listdir(folder) if f.endswith(".yaml")]
    yamls = (yaml.load(open(os.path.join(folder, f))) for f in files)
    pos_files = (f for y, f in zip(yamls, files) if y["cancer"] is label)
    return (os.path.join(folder, f.replace("_metadata.yaml", ".npy.gz")) for f in pos_files)


def get_kmers(ra, k):
    encoded = RaggedArray(bionumpy.encodings.AminoAcidEncoding.from_bytes(
        ra.ravel()), ra.shape)
    kmers = np.lib.stride_tricks.sliding_window_view(encoded.ravel(), k)
    hashes = bionumpy.sequence.kmers.KmerEncoder(k, 20).from_bytes(kmers)
    ra = RaggedArray(hashes, encoded.shape)
    return ra[:, :-(k-1)]


def train_em(kmers, k):
    n_kmers = 20**k
    weights = np.ones(2*n_kmers)+np.random.rand(2*n_kmers)/100
    weights = weights.reshape(2, -1)
    log_ps = np.log(weights/weights.sum(axis=-1, keepdims=True))
    print(log_ps.shape)
    dists = [EM.MultiCategoricalRegularized(lp) for lp in log_ps]
    em = EM.SparseEM(dists, np.log([0.1, 0.9]))
    em.fit(kmers)
    return em


def train_mil(pos_kmers, neg_kmers, k):
    n_kmers = 20**k
    weights = np.ones(2*n_kmers)+np.random.rand(2*n_kmers)/100
    weights = weights.reshape(2, -1)
    log_ps = np.log(weights/weights.sum(axis=-1, keepdims=True))
    print(log_ps.shape)
    dists = [EM.MultiCategoricalRegularized(lp) for lp in log_ps]
    mil = EM.SparseMIL(dists, np.log([0.99, 0.01]))
    mil.fit(pos_kmers, neg_kmers)
    return mil


def get_logodds(counts):
    N = counts.sum()
    return np.log(counts)-np.log(N-counts)


def simple_logistic(pos_kmers, neg_kmers, k):
    pos_counts = np.bincount(pos_kmers.ravel(), minlength=20**k)
    neg_counts = np.bincount(neg_kmers.ravel(), minlength=20**k)
    eta_diff = get_logodds(pos_counts)-get_logodds(neg_counts)
    return np.argpartition(eta_diff, -64)[-64:]


if True or __name__ == "__main__":
    if False:
        folder = "/home/knut/Data/datasets_for_knut/n_kmers_64_repsize_100k_witnessrate_0.0005/immuneml/repertoires/"
        all_kmers = {}
        for label in (True, False):
            seqs = np.concatenate([read_repertoire(filename, filename.replace(".npy.gz", "_metadata.yaml"))
                                   for filename in get_labelled_repertoirs_filenames(folder, label)])
            kmers = get_kmers(seqs, 4)
            kmers.save(f"kmers_{label}-k4.npz")
            all_kmers[label] = kmers
            
        #seqs = {label:
        #        np.concatenate([read_repertoire(filename, filename.replace(".npy.gz", "_metadata.yaml"))
        #                        for filename in get_labelled_repertoirs_filenames(folder, label)])
        #        for label in [True, False]}

        #pos_files = get_positive_repertoires_filenames(folder)
        #all_seqs = np.concatenate([read_repertoire(filename, filename.replace(".npy.gz", "_metadata.yaml")) 
        #                           for filename in pos_files])
        # kmers = {label: get_kmers(seqs[label], 4) for label in seqs}
    all_kmers = {label: RaggedArray.load(f"kmers_{label}-k4.npz") for label in (True, False)}
    # signal = simple_logistic(all_kmers[True], all_kmers[False], 4)
    # model = train_em(all_kmers[True], 4)
    model = train_mil(all_kmers[True], all_kmers[False], 4)
    
