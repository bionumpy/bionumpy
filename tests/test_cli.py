import inspect

import bionumpy as bnp
from bionumpy.cli import CliWrapper
from .util import get_file_name


def mock_function(reads: bnp.datatypes.SequenceEntry) -> bnp.datatypes.SequenceEntry:
    return reads[reads.sequence[:, 0] == 'A']


def test_cli_wrapper(data_path, tmp_path):
    cli_function = CliWrapper()(mock_function)
    output_filename = tmp_path / 'tmp.fq.gz'
    input_filename = data_path / 'big.fq.gz'
    cli_function(input_filename, output=output_filename)
    assert bnp.count_entries(output_filename) < bnp.count_entries(input_filename) // 2


def test_cli_wrapper_annotations():
    cli_function = CliWrapper()(mock_function)
    argspec = inspect.getfullargspec(cli_function)
    print(argspec)
    assert argspec.annotations['reads'] == str
