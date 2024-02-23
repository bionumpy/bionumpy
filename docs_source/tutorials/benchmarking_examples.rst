.. benchmarking_examples

Benchmarking Examples
---------------------

Following is the code for the examples used in the benchmarking. The scripts will usually work both as a command line script,
but also includes a test showing how to use it in python.

BAM Filtering
-------------
Read a BAM file and filter out reads with MAPQ < 60:

.. code-block:: python

        with bnp.open(output[0], 'w') as f:
            for chunk in bnp.open(input[0]).read_chunks():
                f.write(chunk[chunk.mapq == 60])


Unique Intersection
-------------------
Filter one bed file on whether the intervals overlap with another bed file:

.. literalinclude:: /../scripts/unique_intersect_example.py
    :language: python

Jaccard Index
-------------
Calculate the Jaccard index between two or more bed files:

.. literalinclude:: /../scripts/jaccard_all_vs_all_example.py
    :language: python