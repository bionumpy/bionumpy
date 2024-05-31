import numpy as np
from npstructures import RunLengthArray
from .genomic_data import GenomicArray, GenomicIntervals
from .genomic_data.binned_genome import BinnedGenome
from .io.matrix_dump import Matrix
from .sequence.count_encoded import EncodedCounts
from .encoded_array import EncodedRaggedArray
import dataclasses


@dataclasses.dataclass
class Vector:
    data: np.ndarray
    names: list


class Plotter:
    def __init__(self, plt=None):
        self._plt = plt
        self._show = True
        self._tried = False

    @property
    def plt(self):
        if self._plt is None and not self._tried:
            try:
                import matplotlib.pyplot as _plt
                _plt.style.use('Solarize_Light2')
                self._plt = _plt
            except:
                pass
            self._tried = True
        return self._plt


    def set_config(self, **kwargs):
        accepted_configs = {'show'}
        for key, value in kwargs.items():
            assert key in accepted_configs
            if key == 'show':
                self._show = value

    def show(self, f=None):
        if not self._show:
            return
        plt.show()
        # if f is None:
        #     self.plt.show()
        # else:
        #     print(';;;;;;;;;;;;;;;;')
        #     f.show()

    def _conversion(self, data):
        if isinstance(data, GenomicIntervals):
            return data.get_pileup()
        elif isinstance(data, EncodedCounts):
            if len(data.counts.shape) == 2:
                return Matrix(data.counts, col_names=data.alphabet)
            else:
                return Vector(data.counts, names=data.alphabet)
        return data

    def _plot_bars(self, vector):
        fig, ax = self.plt.subplots()
        ax.bar([str(c) for c in vector.names], vector.data)
        self.show()

    def _plot_heatmap(self, matrix):

        fig, ax = self.plt.subplots()
        data = np.asanyarray(matrix.data)
        n_rows, n_cols = data.shape
        im = ax.imshow(data)

        # Show all ticks and label them with the respective list entries

        ax.set_xticks(np.arange(n_cols))
        print(matrix.col_names, type(matrix.col_names))
        if matrix.col_names is not None:
            if isinstance(matrix.col_names, EncodedRaggedArray):
                ax.set_xticklabels(matrix.col_names.tolist())
            else:
                ax.set_xticklabels(matrix.col_names)
        ax.set_yticks(np.arange(n_rows))
        if matrix.row_names is not None:
            if isinstance(matrix.row_names, EncodedRaggedArray):
                ax.set_yticklabels(matrix.row_names.tolist())
            else:
                ax.set_yticklabels(matrix.row_names)

        # Rotate the tick labels and set their alignment.
        self.plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
                      rotation_mode="anchor")

        # Loop over data dimensions and create text annotations.
        if n_rows * n_cols < 100:
            for i in range(n_rows):
                for j in range(n_cols):
                    text = ax.text(j, i, data[i, j],
                                   ha="center", va="center", color="w")
        fig.tight_layout()
        self.show()

    def _plot_list(self, data_list):
        data_list = [self._conversion(data) for data in data_list]
        if isinstance(data_list[0], GenomicArray):
            gc = data_list[0].genome_context
            assert all(isinstance(data, GenomicArray) for data in data_list)
            f, axes = self._subplots_for_genome(gc)
            # f, axes = self.plt.subplots(1, len(gc.chrom_sizes), sharey=True, sharex=True)
            for i, chromosome in enumerate(gc.chrom_sizes.keys()):
                for data in data_list:
                    axes[i].plot(np.asarray(data[chromosome]))
                axes[i].title.set_text(chromosome)
            self.show(f)
        else:
            raise NotImplementedError(f'{type(data)} not supported')

    def _plot_dict(self, data_dict):
        data_dict = {key: self._conversion(data) for key, data in data_dict.items()}
        data = list(data_dict.values())[0]
        if isinstance(data, GenomicArray):
            gc = data.genome_context
            assert all(isinstance(data, GenomicArray) for data in data_dict.values())
            f, axes = self._subplots_for_genome(gc)
            for i, chromosome in enumerate(gc.chrom_sizes.keys()):
                for key, data in data_dict.items():
                    if i == 0:
                        label = key
                    else:
                        label = None
                    axes[i].plot(np.asarray(data[chromosome]), label=label)
                axes[i].title.set_text(chromosome)
            f.legend()
            self.show(f)
        elif isinstance(data, RunLengthArray):
            f, axes = self.plt.subplots()
            for key, data in data_dict.items():
                self._plot_single(data, ax=axes, label=key)
            f.legend()
            self.show()
        else:
            raise NotImplementedError(f'{type(data)} not supported')
        pass

    def _subplots_for_genome(self, gc):
        return self.plt.subplots(len(gc.chrom_sizes), 1, sharey=True, sharex=True)

    def _plot_single(self, data, ax=None, label=None):
        data = self._conversion(data)
        if isinstance(data, (GenomicArray, BinnedGenome)):
            f, axes = self._subplots_for_genome(data.genome_context)
            # f, axes = self.plt.subplots(1, len(data.genome_context.chrom_sizes), sharey=True, sharex=True)
            for i, chromosome in enumerate(data.genome_context.chrom_sizes.keys()):
                axes[i].plot(np.asarray(data[chromosome]))
                axes[i].title.set_text(chromosome)
            self.show(f)
        elif isinstance(data, RunLengthArray):
            if ax is None:
                f, ax = self.plt.subplots()  # self._subplots_for_genome(data.genome_context)
                ax.plot(np.asarray(data), label=label)
                self.show()
            else:
                ax.plot(np.asarray(data), label=label)
        elif isinstance(data, Matrix):
            return self._plot_heatmap(data)
        elif isinstance(data, Vector):
            return self._plot_bars(data)
        elif isinstance(data, np.ndarray):
            f, ax = self.plt.subplots()
            ax.plot(data)
            self.show()
        else:
            raise NotImplementedError(f'{type(data)} not supported')

    def __call__(self, data):
        return self.plot(data)

    def plot(self, data):
        if isinstance(data, list):
            return self._plot_list(data)
        elif isinstance(data, dict):
            return self._plot_dict(data)
        return self._plot_single(data)


def convert_args_to_array(func):
    def convert(a):
        if isinstance(a, list):
            return [np.asarray(e) for e in a]
        elif isinstance(a, dict):
            return {key: np.asarray(val) for key, val in a.items()}
        else:
            return np.asarray(a)

    def new_func(*args, **kwargs):
        args = [convert(a) for a in args]
        kwargs = {key: convert(val) for key, val in kwargs.items()}
        return func(*args, **kwargs)

    return new_func


class PXWrapper:
    def __init__(self, px=None):
        self._px = px
        self._tried = False

    @property
    def px(self):
        if self._px is None and not self._tried:
            try:
                import plotly.express as _px
                self._px = _px
            except:
                pass
            self._tried = True
        return self._px

    def __getattr__(self, name):
        func = getattr(self.px, name)
        return convert_args_to_array(func)

# px = PXWrapper()
# try:
#     import plotly.express as _px
#     px = PXWrapper(_px)
# except:
#     px = None

# try:
#     import matplotlib.pyplot as plt
#     plt.style.use('Solarize_Light2')
plot = Plotter()
# except:
#     plot = None
