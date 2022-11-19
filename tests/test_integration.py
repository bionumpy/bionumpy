import dataclasses
import os
from pathlib import Path

from npstructures import RaggedArray
from npstructures.testing import assert_raggedarray_equal
from bionumpy import AminoAcidEncoding, DNAEncoding, as_encoded_array
from bionumpy.io.files import bnp_open
from bionumpy.io.delimited_buffers import get_bufferclass_for_datatype
from bionumpy.sequence.string_matcher import RegexMatcher
from bionumpy.bnpdataclass import bnpdataclass


def test_csv_import_and_regex_matching():
    seqs_path = Path("seqs.tsv")

    @bnpdataclass
    class SeqAsTSV:
        sequence_aa: AminoAcidEncoding
        sequence: DNAEncoding
        v_call: str
        j_call: str

    buffer_class = get_bufferclass_for_datatype(SeqAsTSV, delimiter='\t', has_header=True)

    sequences_to_write = SeqAsTSV(sequence_aa=as_encoded_array(["AAACCC", "EEAAF"], AminoAcidEncoding),
                                  sequence=as_encoded_array(["AT", "CCA"], DNAEncoding),
                                  v_call=as_encoded_array(['V1-1', 'V2']), j_call=as_encoded_array(['J2', 'J3-2']))

    with bnp_open(seqs_path, buffer_type=buffer_class, mode='w') as file:
        file.write(sequences_to_write)

    file = bnp_open(str(seqs_path), buffer_type=buffer_class)
    sequences = file.read()
    file.close()
    assert isinstance(sequences.sequence_aa, RaggedArray), sequences.sequence_aa

    matcher = RegexMatcher('AA', encoding=AminoAcidEncoding)
    matches = matcher.rolling_window(sequences.sequence_aa, mode='same')
    assert_raggedarray_equal(matches, [[True, True, False, False, False, False],
                                       [False, False, True, False, False]])

    for field in dataclasses.fields(sequences):
        assert_raggedarray_equal(getattr(sequences_to_write, field.name), getattr(sequences, field.name))

    if seqs_path.is_file():
        os.remove(seqs_path)
