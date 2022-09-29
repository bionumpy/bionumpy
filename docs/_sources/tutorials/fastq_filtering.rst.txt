FASTQ filtering
-----------------
Before following this tutorial, we assume you have already followed the introduction part of reading files (see :ref:`reading_files`).

The following in an example of a small script that filters FASTQ reads. This example illustrates the use of multiple functions decorated with `@streamable()`. Each function is designed so that it initially works on one chunk, but with the streamable descorator, we can send chunks from a file and BioNumPy handles the rest for us.

This example also illustrates how to chain multiple functions.

.. code-block:: python

    import numpy as np
    import bionumpy as bnp
    from bionumpy.npdataclassstream import streamable

    @streamable()
    def filter_reads_on_mean_base_quality(reads, minimum_base_quality=20):
        mask = np.mean(reads.quality, axis=-1) > minimum_base_quality
        return reads[mask]

    @streamable()
    def filter_reads_on_minimum_base_quality(reads, min_base_quality=5):
        mask = np.min(reads.quality, axis=-1) > min_base_quality
        return reads[mask]

    def main():
        reads = bnp.open("example_data/big.fq.gz").read_chunks()
        reads = filter_reads_on_mean_base_quality(reads, 10)
        reads = filter_reads_on_minimum_base_quality(reads, 1)

        print("Number of reads after filtering: ", sum(len(r) for r in reads))

