import gzip
import numpy as np
import yaml
from importlib import reload
import bionumpy.kmers
from npstructures import RaggedArray
import bionumpy.encodings
reload(bionumpy.encodings)
reload(bionumpy.kmers)


def read_repertoire(data_filename, meta_filename):
    y = yaml.load(open(meta_filename))
    i = y["field_list"].index("sequence_aas")
    d = np.load(gzip.open(data_filename))
    seqs = [row[i] for row in d]
    return RaggedArray([[ord(c) for c in seq] for seq in seqs])


def get_kmers(ra):
    encoded = RaggedArray(bionumpy.encodings.AminoAcidEncoding.from_bytes(
        seqs.ravel()), seqs.shape)

    kmers = np.lib.stride_tricks.sliding_window_view(encoded.ravel(), 3)
    hashes = bionumpy.kmers.KmerHash(20, 3).hash(kmers)
    ra = RaggedArray(hashes, encoded.shape)
    return ra[:, :-(3-1)]


if True or __name__ == "__main__":
    default_name = "/home/knut/Downloads/datasets_for_knut/n_kmers_64_repsize_100k_witnessrate_0.0005/immuneml/repertoires/7cd3afc8601e4a889da231c6eff7179b.npy.gz"
    filename = default_name
    print(filename)
    seqs = read_repertoire(filename, filename.replace(".npy.gz", "_metadata.yaml"))
    kmers = get_kmers(seqs)
    print(kmers)
