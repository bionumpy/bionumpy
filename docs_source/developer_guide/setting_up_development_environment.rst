
.. _development_setup:

Setting up your development environment
========================================

Follow the steps below to setup an environment that lets you do development on BioNumPy (including making examples or writing documentation).

Step 1 (Optional): Install development version of npstructures
****************************************************************
Development on BioNumPy is often done in parallel with development on NpStructures, meaning that


Step 2: Install BioNumPy locally (not from PyPi)
**************************************************

In order to develop on BioNumPy you will need to clone this repository and install BioNumPy:

.. code-block:: bash

    git clone git@github.com:bionumpy/bionumpy.git

We suggest installing a local editable copy of BioNumPy so that you get BioNumPy in your Python path. It can be a good idea to create a Python virtual environment before doing this. If you installed NpStructures in the previous step, make sure you use the same virtual environment now:

.. code-block:: bash

    cd bionumpy
    pip install -e .


Step 3: Install development dependencies
******************************************
The above installation will install all required dependencies that BioNumPy uses, but for development there are some more dependencies that you will need (e.g. for testing). Install these by running:

.. code-block:: bash

    pip install -r requirements_dev.txt


Step 4: Test your installation
********************************
BioNumPy has several different sets of test (unit tests, property testing, testing of examples and doctesting). All should pass if you've setup everything correctly. You can run all tests by running the `test`-script:

.. code-block:: bash

    ./test

If the above runs without any errors, you are ready to start doing development. Continue by reading the







