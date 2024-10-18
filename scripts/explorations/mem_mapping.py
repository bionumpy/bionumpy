import numpy as np
import bionumpy as bnp
if __name__ == '__main__':
    # Create a new memory-mapped array
    shape = (10)
    dtype = np.uint8
    lengths = np.memmap('myfile.dat',
                        dtype=np.int32, mode='w+', shape=4)
    lengths[:] = [1, 2, 3, 4]

    # Create a new memory-mapped file on disk
    data = np.memmap('data.dat',
                     dtype=dtype, mode='w+', shape=shape)
    data[:] = [ord(c) for c in 'abcdefghij']
    encoded_array = bnp.EncodedArray(
        data,
        bnp.BaseEncoding)  # type: ignore

    print(encoded_array)

    encoded_ragged_array = bnp.EncodedRaggedArray(
        encoded_array, lengths)  # type: ignore

    # Write data to the memory-mapped file
    print(encoded_ragged_array)

    # Flush changes to disk
    data.flush()

    # Clean up when done (the file will remain on disk)
    del data
