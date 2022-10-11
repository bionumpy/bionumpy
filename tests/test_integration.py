import os
from pathlib import Path

import numpy as np
from npstructures import npdataclass, RaggedArray

import bionumpy.encodings
from bionumpy import AminoAcidArray, DNAArray
from bionumpy.files import bnp_open
from bionumpy.delimited_buffers import get_bufferclass_for_datatype
from bionumpy.string_matcher import RegexMatcher


def test_csv_import_and_regex_matching():

    seqs_path = Path("seqs.tsv")

    with open(seqs_path, 'w') as file:
        file.write("""sequence_aa	sequence	v_call	j_call
AAACCC	A	V1-1	J2
EEAAF	CCA	V2	J3-2
""")

    @npdataclass
    class SeqAsTSV:
        sequence_aa: AminoAcidArray
        sequence: DNAArray
        v_call: str
        j_call: str

    sequences = bnp_open(str(seqs_path), buffer_type=get_bufferclass_for_datatype(SeqAsTSV, delimiter='\t'), has_header=True).read()

    assert isinstance(sequences.sequence_aa, RaggedArray)

    matcher = RegexMatcher('AA', encoding=bionumpy.encodings.AminoAcidEncoding)
    matches = matcher.rolling_window(sequences.sequence_aa, mode='same')

    assert np.array_equal(matches, [[True, True, False, False, False, False], [False, False, True, False, False]])

    if seqs_path.is_file():
        os.remove(seqs_path)
