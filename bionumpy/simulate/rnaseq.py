from itertools import chain

def get_transcript_copies(sequences, sequence_counts):
    indices = list(chain(*[[i] * count for i, count in enumerate(sequence_counts)]))
    return sequences[indices]

def fragment_transcript_copies(sequences, fragment_size):
    fragments = []
    for sequence in sequences:
        for i in range(0, len(sequence) - fragment_size + 1, fragment_size):
            fragments.append(sequence[i:i + fragment_size])
    return fragments

