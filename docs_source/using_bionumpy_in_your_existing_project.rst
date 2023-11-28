.. _using_bionumpy_in_your_existing_project:

Using BioNumPy in your existing project
--------------------------------------

If you already have a Python project that reads biological datasets (e.g. fastq, bam or vcf files), you can often easily add BioNumPy to your project to speed up reading and parsing of such files.

All you need is to install BioNumPy, import it, and read your file with `bnp.open`:

.. testcode::

    import numpy as np
    import bionumpy as bnp

    # open your file, bnp.open automatically detects the file format
    f = bnp.open("example_data/variants.vcf")
	# a chunk is an efficient representation of a chunk of many lines
    for chunk in f.read_chunks():
		# we can iterate over the entries
        for single_entry in chunk.to_iter():
			print(single_entry)
			# and we can access things like, chromosome, position and so on
			position = single_entry.position


The above example shows how you can read any file and iterate over the lines. This way, one can in most cases easily use BioNumPy as an inplace replacement for the existing code that iterates over lines from biological data files. However, iterating over lines is not very efficient, so we can do better by working directly on chunks of the file. See :ref:`reading_files` for more information.

Using BioNumPy with Pandas
============================

Another way that BioNumPy can be used directly in existing projects is if you are using Pandas.

It is straightforward to convert a chunk of a file or a whole file to a Pandas DataFrame. This time, we read a fastq file:

.. testcode::

	import bionumpy as bnp
	f = bnp.open("example_data/big.fq.gz")
	for chunk in f.read_chunks():
		df = chunk.topandas()
		print(df)

.. testoutput::

	asef

Continue to see an overview of :ref:`what you can do with bionumpy<what_can_you_do>`.

