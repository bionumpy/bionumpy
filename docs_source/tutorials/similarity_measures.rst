.. _similarity_measures_tutorial:

Computing the similarity between two BED-files
-----------------------------------------------

BioNumPy supports computing the Forbes and the Jaccard similarity measures between two bed files.
The bed-files needs to be sorted, but BioNumPy can do the sorting.

In the following example, we read two bed files (transcription factor binding sites for CTCF and ZNF263) and compute the similarity:

.. literalinclude:: /../scripts/forbes_example.py



