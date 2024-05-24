from bionumpy.cli import CliWrapper
from bionumpy.datatypes import SequenceEntry


def mock_function(reads: SequenceEntry, letter: str='A') -> SequenceEntry:
    '''Subsample a set of reads'''

    return reads[reads.sequence[:, 0] == letter]


if __name__ == "__main__":
    cli_function = CliWrapper()(mock_function)
    import typer

    typer.run(cli_function)
