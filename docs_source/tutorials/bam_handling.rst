BAM Handling
-----------------
Before following this tutorial, we assume you have already followed the introduction part of reading files (see :ref:`reading_files`).

The following example shows how to handle bam files. It checks whether base pairs pairs that are softclipped in alignments differ from the base pairs in the rest of the reads. It reads in a .bam file and uses the cigar_op attribute to find softclipped reads. 


.. code-block:: python

    import bionumpy as bnp
    from npstructures import ragged_slice
    
    # Open the aligments file
    alignments = bnp.open("example_data/test.bam").read()
    
    # Extract the first cigar operation for each alignment
    start_cigar = alignments.cigar_op[..., 0]
    
    # Get aligments that start with soft-clip
    start_clipped_alignments = alignments[start_cigar == "s"]
    
    # Get the number of softclipped
    n_clipped_bases = start_clipped_alignments.cigar_length[..., 0]
    
    # Extract clipped bases
    clipped_bases = ragged_slice(start_clipped_alignments.sequence,
                                 ends=n_clipped_bases)
    
    # Count bases in softclipped regions
    print(bnp.count_encoded(clipped_bases.ravel()))
    
    # Count bases in whole reads
    print(bnp.count_encoded(alignments.sequence.ravel()))

