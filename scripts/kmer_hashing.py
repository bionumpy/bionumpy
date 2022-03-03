import sys
from bionumpy.parser import BufferedNumpyParser
from bionumpy.kmers import TwoBitHash

parser = BufferedNumpyParser.from_filename(sys.argv[1], chunk_size=10000000)
hasher = TwoBitHash(k=31)
for chunk in parser.get_chunks():
    #Can move data to gpu already here
    #chunk is FileBuffer object

    sequences = chunk.get_sequences()
    kmers = hasher.get_kmer_hashes(sequences)
    # kmers can be passed to counter
