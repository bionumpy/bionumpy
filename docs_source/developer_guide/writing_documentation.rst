.. _writing_documentation:

Writing documentation
=========================

To write documentation, edit or create new files in the `docs_source` directory.

To build and test the documentation locally, run:

.. code-block:: bash

    cd docs_source
    make html

This will create html files in the `docs_source/_build` directory. Open the index.html file in a web browser to see the build.


Note: If you make big changes to the documentation, you will need to clean old html files:

.. code-block:: bash

    make clean


Structure of documentation
------------------------------
* `source`: Put documentation of BioNumPy concepts in this directory. You will need to manually add a link from `index.rst` to new files created.
* `tutorials`: Put tutorials here. These will automatically be indexed.
* `developer_guide`: Put developer guides here. You will need to manually add a link from `index.rst` to new files created.


Writing code in the documentation
----------------------------------
Code written in the documentation will be tested automatically with doctest. For instance, if writing the following, doctest will check whether the output from the Python code matches the provided output:

    >>> 2 + 2
    4

You can also write code that is tested by using the testcode/testoutput directives. `Read this guide to see how it is done <https://www.sphinx-doc.org/en/master/usage/extensions/doctest.html#directive-testcode>`_

A thid option is to include code from e.g. the scripts directory directly into the documentation:

.. code-block::

    .. literalinclude:: /../scripts/your_example.py


Check that all code you write passes the tests by running:

.. code-block:: bash
    make doctest

You may get a lot of warnings, but the output should end with `Build succeeded`.
