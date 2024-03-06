import numpy as np
from npstructures import RaggedArray

from .. import SequenceEntryWithQuality
from ..genomic_data import GenomicSequence
from ..simulate import simulate_sequences
import pytest
import bionumpy as bnp
from ..simulate.intervals import simulate_intervals


@pytest.fixture
def rng():
    return np.random.default_rng(0)


@pytest.fixture
def chromosome_sizes():
    return {"chr1": 100, "chr2": 200, "chr15": 300}


@pytest.fixture
def genome_sequence(rng, chromosome_sizes):
    return simulate_sequences("acgt", chromosome_sizes, rng=rng)


@pytest.fixture
def intervals(chromosome_sizes):
    return simulate_intervals(chromosome_sizes, 10, 20)


@pytest.fixture
def intervals2(chromosome_sizes):
    return simulate_intervals(chromosome_sizes, 15, 30)


@pytest.fixture
def reads(genome_sequence, intervals):
    gs = GenomicSequence.from_dict({entry.name: entry.sequence for entry in genome_sequence.tolist()})
    seqs = gs[intervals]
    return SequenceEntryWithQuality([f'seq{i}' for i in range(len(seqs))],
                                    seqs, RaggedArray(np.full(sum(seqs.lengths), 40), seqs.lengths))


@pytest.fixture
def pileup(genome: bnp.Genome, intervals):
    return genome.get_intervals(intervals).get_pileup()


@pytest.fixture()
def genome(chromosome_sizes):
    return bnp.Genome.from_dict(chromosome_sizes)


@pytest.fixture
def fasta_file(tmp_path, genome_sequence):
    fasta_file = tmp_path / "genome.fa"
    with bnp.open(fasta_file, "w") as f:
        f.write(genome_sequence)
    return fasta_file


@pytest.fixture
def bed_file(tmp_path, intervals):
    bed_file = tmp_path / "intervals.bed"
    with bnp.open(bed_file, "w") as f:
        f.write(intervals)
    return bed_file


@pytest.fixture
def bed_file2(tmp_path, intervals):
    bed_file = tmp_path / "intervals.bed"
    with bnp.open(bed_file, "w") as f:
        f.write(intervals)
    return bed_file


@pytest.fixture
def fastq_file(tmp_path, reads):
    fastq_file = tmp_path / "reads.fastq"
    with bnp.open(fastq_file, "w") as f:
        f.write(reads)
    return fastq_file


@pytest.fixture
def bedgraph_file(tmp_path, pileup):
    pileup = pileup.get_data()
    bedgraph_file = tmp_path / "pileup.bdg"
    with bnp.open(bedgraph_file, "w") as f:
        f.write(pileup)
    return bedgraph_file
