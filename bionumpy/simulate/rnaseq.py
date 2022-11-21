from itertools import chain

def get_transcript_copies(sequences, sequence_counts):
    indices = list(chain(*[[i] * count for i, count in enumerate(sequence_counts)]))
    return sequences[indices]