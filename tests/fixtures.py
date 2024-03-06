from bionumpy import Bed6
from bionumpy.datatypes import VCFEntry, SequenceEntryWithQuality

data = {
    "bed": Bed6.from_entry_tuples([
        ("chr1", 1, 3, ".", 0, "-"),
        ("chr1", 40, 60, ".", 1, "+"),
        ("chr20", 400, 600, ".", 2, "+")]),
    "vcf2": VCFEntry.from_entry_tuples([
        ("chr1", 88361, "rs4970378", "A", "G", ".", ".", "."),
        ("chr1", 887559, "rs3748595", "A", "CAA", ".", ".", "."),
        ("chr2", 8877, "rs3828047", "AGG", "C", ".", ".", ".")]),
    "fastq": SequenceEntryWithQuality.from_entry_tuples([
        ("headerishere", "CTTGTTGA", "".join("!" for _ in "CTTGTTGA")),
        ("anotherheader", "CGG", "".join("~" for _ in "CGG"))]), }
'''
    "vcf": [
        VCFEntry.single_entry("chr1",	88361, "rs4970378",	"A",	"G", ".", ".", "."),
        VCFEntry.single_entry("chr1",	887559, "rs3748595",	"A",	"C", ".", ".", "."),
        VCFEntry.single_entry("chr2",	8877, "rs3828047",	"A",	"G", ".", ".", ".")],

    "fasta": [
        SequenceEntry.single_entry("header", "CTTGTTGA"),
        SequenceEntry.single_entry("header2", "CGG")],
    "multiline_fasta": [
        SequenceEntry.single_entry("header", "CTTGCCGCCTCC"),
        SequenceEntry.single_entry("header2", "CCCCCCGGGCCCTTT")],
    "gfa_sequence": [
        SequenceEntry.single_entry("id1", "AACCTTGG"),
        SequenceEntry.single_entry("id4", "ACTG")],
    "sam": [
        SAMEntry.single_entry("SRR1524970.144283", 16, "test_ref", 1705, 255, "25M",	"*", 0, 0, "TGCTGATGAAGCAGAACAACTTTAA", "]YG[^baaaa^W`ab]]````aaba", ''),
        SAMEntry.single_entry("SRR1524970.316478", 16, "test_ref", 1705, 255, "24M", "*", 0, 0, "TGCTGATGAAGCAGAACAACTTTA", 	"`\X_`aaaaaY]``b_aa_aaaaa", '')],
    "bed12": [
        Bed12.single_entry("chr21", 10079666,  10120808,   "uc002yiv.1", 0, "-", 10081686, 10120608,  "0",     4,   [[528,91,101,215]], [[0,1930,39750,40927]]),
        Bed12.single_entry("chr21", 10080031,  10081687,   "uc002yiw.1",  0,  "-",  10080031,  10080031,  "0",     2,   [[200,91]],    [[0,1565]])],
    'wig': BedGraph.from_entry_tuples([
        ('chr1',	0, 9800,	-0),
        ('chr1',	9800,	9871,	0.36612),
        ('chr1',	9871,	9872,	0.17042)])
}'''


