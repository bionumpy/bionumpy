========
BioNumPy
========

.. image:: https://img.shields.io/pypi/v/bionumpy.svg
        :target: https://pypi.python.org/pypi/bionumpy

.. image:: https://github.com/knutdrand/bionumpy/actions/workflows/python-install-and-test.yml/badge.svg
        :target: https://github.com/bionumpy/bionumpy/actions/
        :alt: Build and test status

Documentation: `https://bionumpy.github.io/bionumpy/ <https://bionumpy.github.io/bionumpy/>`_


What is BioNumPy?
-----------------
BioNumPy is a toolkit, built on top of NumPy, for enabling array programming on biological data in Python. BioNumPy aims to make it easy to read common bioinformatics file formats efficiently into NumPy-like data structures that enable efficient operations and analysis of the data. Working in BioNumPy should feel much like working in NumPy.


Why BioNumPy?
-------------
* There are no existing packages in Python for getting biological data sets efficiently into NumPy-like data structures.
* Current packages for working with biological data do not use NumPy in an efficient way (e.g. individual sequences are stored as separate NumPy arrays, not together in shared arrays).


Getting started
----------------

1. Install with pip:

* Free software: MIT license
* Documentation: https://bionumpy.github.io/bionumpy/.

$ pip install bionumpy

2. Check out the tutorials and getting started guide in the `documentation <https://bionumpy.github.io/bionumpy/>`_.


Features
------------

The features of BioNumPy can roughly be divided into two:

1. Reading biological data sets (e.g. fasta, vcf, bed) into NumPy-like objects
2. Analysing, modifying, filtering these NumPy-like objects in meaningful ways

BioNumPy also supports writing most data types to file.



Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
