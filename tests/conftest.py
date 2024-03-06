from pathlib import Path

import pytest


@pytest.fixture()
def data_path():
    return Path(__file__).parent.parent / 'example_data'
