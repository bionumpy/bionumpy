
Encodings
~~~~~~~~~
A central concept in BioNumPy is the *encoding* of data, such as DNA sequence, base qualities, kmers, etc, to ...



The main point of BioNumPy is to leverage the computational power of NumPy for biological data. A key element to this is to represent different types of biological data as numbers (in numpy arrays). The basic mechanism for doing this is by Encoding classes that can encode data types into numbers, and decode them back in to the data type. A driving idea in BioNumPy is to make this happen automatically under the hood, so that a user can choose to ignore the inner workings and (in most cases) relate to sequence data sets in the same way as one would with standard numerical numpy data structures.
