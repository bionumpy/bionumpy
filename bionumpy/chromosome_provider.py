import numpy as np
from itertools import groupby
from .sequences import Sequence
import logging

logger = logging.getLogger(__name__)


class GroupedData:
    pass


class GroupedStream:
    def __init__(self, stream):
        self._stream = stream

    def __iter__(self):
        return iter(self._stream)


class GroupedDict(dict):
    pass


class ChromosomeProvider:
    """Base class for objects that can provide data grouped by chromosomes, 
    either by supplying a stream of `(chromosome_name, data)` tuples, or by
    providing a getitem method so `obj[chromosome_name]` returns the data for the 
    chromsome. 

    This class works together with @ChromosomeMap/@bnp_broadcast so that a function
    decorated with either of those will apply the function to the data for each chromosome
    if a `ChromosomeProvider` instance is given as argument.
    """
    _default_value = None

    @staticmethod
    def get_chrom_name(char_array):
        return "".join(chr(c) for c in char_array).replace("\x00", "")

    @staticmethod
    def _get_chromosome_changes(chromosomes):
        return (
            np.flatnonzero(np.any(chromosomes[1:] != chromosomes[:-1], axis=-1))
            + 1
        )


class ChromosomeStreamProvider(ChromosomeProvider):
    """
    Provides chromosome data as a stream of `(chromosome_name, data)` tuples. 
    Just a wrapper around a stream to make it recognizable to @ChromosomeMap
    """
    def __init__(self, stream):
        self._stream = stream

    def __iter__(self):
        return iter(self._stream)


class ChromosomeFileStreamProvider(ChromosomeStreamProvider):
    """
    Splits the data devided in file buffers into one chunk per chromsome. 
    """

    def __init__(self, file_buffers, default_value=None):
        self._buffers = file_buffers
        self._default_value = default_value

    def __iter__(self):
        chrom_func = lambda b: self.get_chrom_name(b.chromosome[0])
        chrom_grouped = groupby(self._buffers, chrom_func)
        overlay = []
        for start_chrom, group in chrom_grouped:
            group = list(group)
            last_buffer = group[-1]
            chrom_changes = self._get_chromosome_changes(last_buffer.chromosome)
            if (
                len(overlay)
                and self.get_chrom_name(overlay[0][0].chromosome) != start_chrom
            ):
                overlay[0].chromosome = overlay[0].chromosome.array.view(Sequence)
                yield self.get_chrom_name(overlay[0][0].chromosome), overlay[0]
                overlay = []
            if not chrom_changes.size:
                tmp = np.concatenate(overlay + group)
                tmp.chromosome = tmp.chromosome.array.view(Sequence)
                yield (start_chrom, tmp)
                overlay = []
            else:
                l = last_buffer[: chrom_changes[0]]
                tmp = np.concatenate(overlay + group[:-1] + [l])
                tmp.chromosome = tmp.chromosome.array.view(Sequence)
                yield (
                    start_chrom,
                    tmp,
                )  # np.concatenate(overlay + group[:-1] + [last_buffer[:chrom_changes[0]]]))
                overlay = [last_buffer[chrom_changes[-1] :]]
                if len(chrom_changes > 1):
                    chunks = (
                        last_buffer[start:end]
                        for start, end in zip(chrom_changes[:-1], chrom_changes[1:])
                    )
                    for chunk in chunks:
                        chunk.chromosome = chunk.chromosome.array.view(Sequence)
                        yield (self.get_chrom_name(chunk.chromosome[0]), chunk)
            last_chrom = last_buffer[-1].chromosome
        if len(overlay):
            chunk = overlay[0]
            chunk.chromosome = chunk.chromosome.array.view(Sequence)
            yield (self.get_chrom_name(chunk.chromosome[0]), chunk)

#
#         cur_data = []
#         last_chromosome = np.zeros(3, dtype=np.uint8)
#         last_chromosome = None
#         for file_buffer in self._buffers:
#             chromosomes = file_buffer.chromosome
#             if not len(chromosomes):
#                 break
#             if last_chromosome is not None and not self._is_same_chromosome(last_chromosome, chromosomes[0]) :
#                 yield (self.get_chrom_name(last_chromosome), np.concatenate(cur_data))
#                 last_chromosome = chromosomes[0]
#                 cur_data = []
#             data = file_buffer.data
#             chromosome_changes = self._get_chromosome_changes(chromosomes)
#             if len(chromosome_changes)==0:
#                 cur_data.append(data)
#                 continue
#             cur_data.append(data[:chromosome_changes[0]])
#             yield self.get_chrom_name(last_chromosome), np.concatenate(cur_intervals)
#             for start, end in zip(chromosome_changes[:-1], chromosome_changes[1:]):
#                 yield np.get_chrom_name(chromosomes[start]), data[start:end]
#             last_chromosome = chromosomes[-1]
#             cur_data.append(data[chromosome_changes[-1]:])
#         yield self.get_chrom_name(last_chromosome). np.concatenate(cur_data)


class ChromosomeDictProvider(ChromosomeProvider):
    """
    Makes data for each chromosome accesible by item lookup
    """

    def items(self):
        return self._d.items()

    def values(self):
        return self._d.values()

    def items(self):
        return self._d.items()

    def __getitem__(self, key):
        possible_keys = [key]
        if key.startswith("chr"):
            possible_keys.append(key[3:])
        else:
            possible_keys.append("chr" + key)
        for k in possible_keys:
            if k in self._d:
                return self._d[k]
        logger.warning(
            f"Chromosomedict missing data for chrom: {possible_keys}, inserting {self._default_value}"
        )
        return self._default_value


class PureChromosomeDictProvider(ChromosomeDictProvider):
    def __init__(self, *args, **kwargs):
        self._d = dict(*args, **kwargs)


class FullChromosomeDictProvider(ChromosomeDictProvider):
    def __init__(self, buffers, default_value=None):
        self._d = dict(ChromosomeFileStreamProvider(buffers))
        self._default_value = default_value


class LazyChromosomeDictProvider(ChromosomeDictProvider):
    def __init__(self, buffers, default_value=None):
        self._stream = iter(ChromosomeFileStreamProvider(buffers))
        self._default_value = default_value
        self._back_log = {}

    def _find_item(self, key):
        for chromosome, data in self._stream:
            if chromosome == key:
                return data
            else:
                self._back_log[chromosome] = data
        assert False

    def __getitem__(self, key):
        if key in self._back_log:
            ret = self._back_log[key]
            del self._back_log[key]
            return ret
        return self._find_item(key)
