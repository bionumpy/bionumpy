.. _sequences:

==========
Sequences
==========

`bnp.Sequences` objects are numpy arrays structured in a way that allows them to hold many sequences of unequal lenght. Under the hood they are `npstructures.RaggedArray` objects that hold byte-arrays, with an encoding that specifies which characters each byte represents. Most of the time, it is not necessary to think about these inner workings, but one can think of them as lists of strings (with the possibility of performing numpy functions on them). The most common way to get `Sequences` objects is to read a file, but they can also be created from lists of strings using the `bnp.as_sequence_array` function::

    >>> bnp.as_sequence_array(["acgttgta", "gcttca", "gttattc"])
    Sequences(acgttgta, gcttca, gttattc)

