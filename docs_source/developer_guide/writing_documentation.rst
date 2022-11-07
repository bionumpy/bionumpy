.. _writing_documentation:

Writing documentation
=========================

To write documentation, edit or create new files in the docs_source directory.

To build and test the documentaiton locally, run:

.. code-block:: bash

    cd docs_source
    make html

This will create html files in the `docs_source/_build` directory. Open the index.html file in a web browser to see the build.


Note: If you make big changes to the documentation, you will need to clean old html files:

.. code-block:: bash

    make clean

