import bionumpy as bnp
from bionumpy.plotting import Plotter
from .genomic_fixtures import *
import numpy as np
import pytest

try:
    import matplotlib.pyplot as plt
except:
    plt = None

skipplt = pytest.mark.skipif(plt is None, reason='skippingmatplotlib tests')


@pytest.fixture
def matrix():
    return bnp.Matrix(np.arange(6).reshape(2, 3), ['hei' , 'ja'], ['a', 'b', 'c'])


@skipplt
def test_plotgenomic_array(track):
    plotter = Plotter(plt)
    plotter.set_config(show=False)
    plotter.plot(track)


@skipplt
def test_plotgenomic_array2(track):
    b = track+2
    plotter = Plotter(plt)
    plotter.set_config(show=False)
    plotter.plot([track, b])


@skipplt
def test_plotgenomic_array3(track):
    b = track+2
    plotter = Plotter(plt)
    plotter.set_config(show=False)
    plotter.plot({'a': track, 'b': b})


@skipplt
def test_matrix(matrix):
    plotter = Plotter(plt)
    plotter.set_config(show=False)
    plotter.plot(matrix)
