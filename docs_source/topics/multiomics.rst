Multiomics
===========
Molecular biology data come in range of modalities - DNA, RNA and protein sequences, as well as a broad range of annotations on such sequences, including DNA variation, DNA function, DNA metyhlation, histone modifications, protein binding to DNA, gene expression and much more. Since the philosohpy of bionumpy is to provide a set of fundamental, broad-purpose data representations and operations (rather than a catalog of standard analyses), bionumpy lends itself very well to the myriad ways in which one might want to combine biosequence modalities to answer varied biological or medical questions.

Example: compute GC content inside genes
----------------------------------------
A basic example is to compute the GC content of DNA sequence annotate as constituting gene regions. This inolves two modalities:

1. DNA sequence of a genome
2. A set of genome intervals annotated as genes.

With bionumpy, DNA sequence of a genome can be directly read in from a fasta file, while regions annotated as genes for the same genome is read in from e.g. a bed or GTF file. The DNA sequence within gene regions is easily extracted by simply indexing the genome sequence bionumpy array by the interval bionumpy array of regions annotated as genes. The following tutorial shows how to do this both at small and large scale:

:ref:`Tutorial: GC content inside genes<gc_content>`

