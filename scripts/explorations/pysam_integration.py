import time
from typing import Iterable, Union

import numpy as np
import pysam
import bionumpy as bnp
from bionumpy import Genome, Interval
from bionumpy.alignments import alignment_to_interval
from bionumpy.arithmetics import sort_intervals
from bionumpy.util.testing import assert_bnpdataclass_equal


def extract_bam_file():
    ...


CIGAR_OPS = "MIDNSHP=X"


def alignments_to_bam_entries(alignments: Iterable[pysam.AlignedSegment], none_on_empty=False) -> Union[
    bnp.BamEntry, None]:
    cols = alignments_to_cols(alignments)
    if not cols:
        return bnp.BamEntry.empty() if not none_on_empty else None
    return bnp.BamEntry(*cols)


def alignments_to_cols(alignments):
    attributes = ((read.reference_name, read.query_name, read.flag, read.reference_start,
                   read.mapping_quality, [CIGAR_OPS[c[0]] for c in read.cigartuples],
                   [c[1] for c in read.cigartuples],
                   read.seq, read.qual) for read in alignments)
    cols = [list(col) for col in zip(*attributes)]
    return cols


class IndexedBamFile:
    def __init__(self, filename):
        self._samfile = pysam.AlignmentFile(filename, 'rb')

    def fetch(self, chromosome, start=None, stop=None, none_on_empty=False):
        return alignments_to_bam_entries(self._samfile.fetch(str(chromosome), start=int(start), stop=int(stop)),
                                         none_on_empty=none_on_empty)

    def _fetch_sorted_intervals(self, intervals):
        cur_chromosome, last_stop, last_start = None, None, None
        all_alignments = []
        for interval in intervals:
            alignments = self.fetch(interval.chromosome, start=interval.start, stop=interval.stop, none_on_empty=True)
            if alignments is None:
                continue
            if interval.chromosome == cur_chromosome:
                assert interval.start >= last_start, f'Intervals must be sorted {(interval.start, last_start)}'
                alignments = alignments[alignments.position >= last_stop]
            all_alignments.append(alignments)
            cur_chromosome, last_stop, last_start = interval.chromosome, interval.stop, interval.start
        return np.concatenate(all_alignments)

    def __getitem__(self, item: Interval):
        '''
        Extract all reads that overlap with any interval in the input item.
        Parameters
        ----------
        item

        Returns
        -------

        '''

        item = sort_intervals(item)
        return self._fetch_sorted_intervals(item)


def multi_intervals(bam_filename='../../example_data/ctcf_chr21-22.bam',
                    interval_filename='../../example_data/ctcf.bed.gz',
                    genome_filename='../../example_data/hg38.chrom.sizes'):
    ib = IndexedBamFile(bam_filename)
    intervals = bnp.open(interval_filename).read()
    t = time.time()
    alignments = ib[intervals]
    print('T =', time.time() - t)
    t = time.time()
    my_entries = old_subset(interval_filename, bam_filename, genome_filename)
    print('T =', time.time() - t)
    assert len(alignments) == len(my_entries)


def main():
    import pooch
    bed_filename = pooch.retrieve('https://www.encodeproject.org/files/ENCFF786YUS/@@download/ENCFF786YUS.bed.gz', None)
    genome_filename = pooch.retrieve('https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes', None)
    bam_filename = '/home/knut/Downloads/ENCFF321VRF.bam'
    multi_intervals(bam_filename, bed_filename, genome_filename)


def old_subset(interval_filename, bam_filename='../../example_data/ctcf_chr21-22.bam',
               genome_filename='../../example_data/hg38.chrom.sizes'):
    genome = Genome.from_file(genome_filename, filter_function=None)
    genome_mask = genome.read_intervals(interval_filename).get_mask()
    chunks = bnp.open(bam_filename).read_chunks()
    all_intervals = []
    for my_entries in chunks:
        mask = genome_mask[alignment_to_interval(my_entries)].any(axis=-1)
        all_intervals.append(my_entries[mask])
    return np.concatenate(all_intervals)


if __name__ == '__main__':
    main()
