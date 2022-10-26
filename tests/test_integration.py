import os
from pathlib import Path

import numpy as np
from npstructures import RaggedArray
from npstructures.testing import assert_raggedarray_equal
import bionumpy.encodings
from bionumpy import AminoAcidEncoding, DNAEncoding
from bionumpy.io.files import bnp_open
from bionumpy.io.delimited_buffers import get_bufferclass_for_datatype
from bionumpy.string_matcher import RegexMatcher
from bionumpy.bnpdataclass import bnpdataclass


def test_csv_import_and_regex_matching():

    seqs_path = Path("seqs.tsv")

    with open(seqs_path, 'w') as file:
        file.write("""\
sequence_aa	sequence	v_call	j_call
AAACCC	A	V1-1	J2
EEAAF	CCA	V2	J3-2
""")

    @bnpdataclass
    class SeqAsTSV:
        sequence_aa: AminoAcidEncoding
        sequence: DNAEncoding
        v_call: str
        j_call: str

    buffer_class = get_bufferclass_for_datatype(SeqAsTSV, delimiter='\t', has_header=True)
    sequences = bnp_open(str(seqs_path), buffer_type=buffer_class).read()
    assert isinstance(sequences.sequence_aa, RaggedArray), sequences.sequence_aa

    matcher = RegexMatcher('AA', encoding=AminoAcidEncoding)
    matches = matcher.rolling_window(sequences.sequence_aa, mode='same')
    assert_raggedarray_equal(matches, [[True, True, False, False, False, False],
                                       [False, False, True, False, False]])
    
    if seqs_path.is_file():
        os.remove(seqs_path)
