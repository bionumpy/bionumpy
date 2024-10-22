import dataclasses
import time
from typing import Iterable, Union, List

import numpy as np
import pysam
from npstructures import RaggedArray

import bionumpy as bnp
from bionumpy import Genome, Interval, BamEntry, EncodedArray, EncodedRaggedArray, BaseEncoding
from bionumpy.alignments import alignment_to_interval
from bionumpy.arithmetics import sort_intervals
from bionumpy.encodings import CigarOpEncoding
from bionumpy.util.testing import assert_bnpdataclass_equal


def extract_bam_file():
    ...


CIGAR_OPS = "MIDNSHP=X"


def alignments_to_bam_entries(alignments: Iterable[pysam.AlignedSegment], none_on_empty=False) -> Union[
    bnp.BamEntry, None]:
    cols = alignments_to_cols(alignments)
    return cols.to_entry()
    if not cols:
        return bnp.BamEntry.empty() if not none_on_empty else None
    return bnp.BamEntry(*cols)


@dataclasses.dataclass
class BamAccumulator:
    reference_name: List[str]
    query_name: List[str]
    flag: List[int]
    reference_start: List[int]
    mapping_quality: List[int]
    cigar_ops_data: List[int]
    cigar_ops_lengths: List[int]
    cigar_len_data: List[int]
    cigar_len_lengths: List[int]
    seq_data: List[str]
    seq_lengths: List[int]
    qual_data: List[int]
    qual_lengths: List[int]

    def add_alignment(self, read: pysam.AlignedSegment):
        accumulator = self
        accumulator.reference_name.append(read.reference_name)
        accumulator.query_name.append(read.query_name)
        accumulator.flag.append(read.flag)
        accumulator.reference_start.append(read.reference_start)
        accumulator.mapping_quality.append(read.mapping_quality)
        cigar_ops, cigar_lengths = zip(*read.cigartuples)
        accumulator.cigar_ops_data.extend(cigar_ops)
        accumulator.cigar_ops_lengths.append(len(cigar_ops))
        accumulator.cigar_len_data.extend(cigar_lengths)
        accumulator.cigar_len_lengths.append(len(cigar_lengths))
        seq = read.seq
        accumulator.seq_data.append(seq)
        accumulator.seq_lengths.append(len(seq))
        qual = read.qual
        accumulator.qual_data.append(qual)
        accumulator.qual_lengths.append(len(qual))

    def to_entry(self):
        if not len(self.reference_name):
            return None
        seq_data = np.frombuffer(bytes(''.join(self.seq_data), 'utf-8'), dtype=np.uint8)
        qual_data = np.frombuffer(bytes(''.join(self.qual_data), 'utf-8'), dtype=np.uint8)
        return BamEntry(
            np.array(self.reference_name),
            np.array(self.query_name),
            np.array(self.flag),
            np.array(self.reference_start),
            np.array(self.mapping_quality),
            EncodedRaggedArray(EncodedArray(
                np.array(self.cigar_ops_data), CigarOpEncoding),
                np.array(self.cigar_ops_lengths)),
            RaggedArray(np.array(self.cigar_len_data), np.array(self.cigar_len_lengths)),
            EncodedRaggedArray(EncodedArray(seq_data, BaseEncoding),
                               np.array(self.seq_lengths)),
            RaggedArray(qual_data, np.array(self.qual_lengths)))



def alignments_to_cols(alignments, min_start=0, accumulator=None):
    accumulator = accumulator or BamAccumulator(*(list() for _ in range(13)))
    for read in alignments:
        if read.reference_start < min_start:
            continue
        accumulator.add_alignment(read)
    return accumulator
    #return accumulator.to_entry()


def alignments_to_cols_old(alignments):
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

    def _fetch_from_sorted_intervals(self, intervals):
        cur_chromosome, last_stop, last_start = None, None, None
        accumulator = BamAccumulator(*(list() for _ in range(13)))
        for interval in intervals:
            if interval.chromosome == cur_chromosome:
                assert interval.start >= last_start, f'Intervals must be sorted {(interval.start, last_start)}'
            else:
                last_stop = 0
            alignments_to_cols(self._samfile.fetch(str(interval.chromosome),
                                                   start=int(interval.start),
                                                   stop=int(interval.stop)),
                               min_start=last_stop, accumulator=accumulator)
            #alignments = self.fetch(interval.chromosome, start=interval.start, stop=interval.stop, none_on_empty=True)
            cur_chromosome, last_stop, last_start = interval.chromosome, interval.stop, interval.start
            continue
            if alignments is None:
                continue

            all_alignments.append(alignments)

        return accumulator.to_entry()
        return np.concatenate(all_alignments)

    def get_all_overlapping(self, intervals: Interval) -> BamEntry:
        item = sort_intervals(intervals)
        return self._fetch_from_sorted_intervals(item)


    def __getitem__(self, item: Interval):
        '''
        Extract all reads that overlap with any interval in the input item.
        Parameters
        ----------
        item

        Returns
        -------

        '''
        return self.get_all_overlapping(item)




def multi_intervals(bam_filename='../../example_data/ctcf_chr21-22.bam',
                    interval_filename='../../example_data/ctcf.bed.gz',
                    genome_filename='../../example_data/hg38.chrom.sizes'):
    ib = IndexedBamFile(bam_filename)
    intervals = bnp.open(interval_filename).read()
    t = time.time()
    alignments = ib[intervals]
    print('T =', time.time() - t)
    return alignments
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

def _test():
    multi_intervals()

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
