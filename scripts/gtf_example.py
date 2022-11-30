import bionumpy as bnp
from bionumpy.io.strops import str_equal, split


def find_closest_gene(annotation_file, narrow_peaks_file):
    gtf_entries = bnp.open(annotation_file).read()
    transcript = gtf_entries.get_transcripts()
    peaks = bnp.open(narrow_peaks_file).read()
    print(transcript.strand)
    pos_transcripts = transcript[transcript.strand.ravel() == "+"]
    print(pos_transcripts)

def test():
    find_closest_gene("example_data/small.gtf", "example_data/peaks.narrowPeak")


if __name__ == "__main__":
    test()
