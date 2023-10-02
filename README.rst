========
BioNumPy
========

.. image:: https://img.shields.io/pypi/v/bionumpy.svg
        :target: https://pypi.python.org/pypi/bionumpy

.. image:: https://github.com/bionumpy/bionumpy/actions/workflows/python-install-and-test.yml/badge.svg
        :target: https://github.com/bionumpy/bionumpy/actions/
        :alt: Build and test status

.. image:: https://github.com/bionumpy/bionumpy/actions/workflows/benchmarking_big.yml/badge.svg
        :target: https://github.com/bionumpy/bionumpy/blob/benchmarks/benchmarks/report_big.md
        :alt: Benchmarks


Documentation: `https://bionumpy.github.io/bionumpy/ <https://bionumpy.github.io/bionumpy/>`_


What is BioNumPy?
-----------------
BioNumPy is a Python library, built on top of NumPy, for enabling array programming on biological datasets in Python.
BioNumPy aims to make it easy to read common bioinformatics file formats efficiently into NumPy-like data structures
that enable efficient operations and analysis of the data. Working in BioNumPy should feel much like working in NumPy.


Getting started
----------------

1. Install with pip:

	$ pip install bionumpy

2. Import BioNumPy and read your data, e.g.:

    >>> import bionumpy as bnp
    >>> import numpy as np
    >>> f = bnp.open("example_data/big.fq.gz")
    >>> # chunk contains the sequences of reads and NumPy-functions can be used
    >>> for chunk in f:
    ...      print(np.sum(chunk.sequence == "G"))
    53686

Check out the getting started guide and various tutorials in the `documentation <https://bionumpy.github.io/bionumpy/>`_.

You can also find the latest version of the `BioNumPy manuscript here <https://github.com/bionumpy/bionumpy/blob/manuscript/bionumpymanuscript.pdf>`_.




