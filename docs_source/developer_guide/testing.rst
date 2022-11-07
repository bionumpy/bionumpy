.. _testing:

Testing
---------

Tests for BioNumPy are automatically run when you push to the dev branch or make a pull request, but we recommend that you also run all tests locally. This can be done by running the test-script:

.. code-block:: bash

    ./test

The output will tell you if anything fails.

The above command will run four different types of tests:

* Unit-tests: Quick to run, will test small units of the code)
* Property-tests: Slow to run, will test properties of the code with big sets of random data)
* Testing of examples: pytest is run on all files with pattern script/*_example.py
* Doctesting: This runs doctest on all docstrings.

When doing development, you may only want to run one of these manually to save some time (e.g. only the doctests if you write documentation, or only uni-tests if you develop some code):

Run only unit-tests:

.. code-block:: bash

    pytest

Run doctests:

.. code-block:: bash

    cd docs_source
    make doctest

Run tests in only a single example (or file):

.. code-block:: bash

    pytest example/our_example.py
