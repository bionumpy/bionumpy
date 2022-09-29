from nimfa.methods.seeding import nndsvd
import numpy as np
import torch
import torch.nn
from torch import nn
import sys


class NMF:
    def __init__(self, V, rank, max_iterations=200000, tolerance=1e-8, test_conv=1000, gpu_id=0, generator=None,
                 init_method='nndsvd', floating_point_precision='single', min_iterations=2000):

        """
        Run non-negative matrix factorisation using GPU. Uses beta-divergence.
        Args:
          V: Matrix to be factorised
          rank: (int) number of latent dimensnions to use in factorisation
          max_iterations: (int) Maximum number of update iterations to use during fitting
          tolerance: tolerance to use in convergence tests. Lower numbers give longer times to convergence
          test_conv: (int) How often to test for convergnce
          gpu_id: (int) Which GPU device to use
          generator: random generator, if None (default) datetime is used
          init_method: how to initialise basis and coefficient matrices, options are:
            - random (will always be the same if set generator != None)
            - NNDSVD
            - NNDSVDa (fill in the zero elements with the average),
            - NNDSVDar (fill in the zero elements with random values in the space [0:average/100]).
          floating_point_precision: (string or type). Can be `double`, `float` or any type/string which
              torch can interpret.
          min_iterations: the minimum number of iterations to execute before termination. Useful when using
              fp32 tensors as convergence can happen too early.
        """
        # torch.cuda.set_device(gpu_id)

        if floating_point_precision == 'single':
            self._tensor_type = torch.FloatTensor
            self._np_dtype = np.float32
        elif floating_point_precision == 'double':
            self._tensor_type = torch.DoubleTensor
            self._np_dtype = np.float64
        else:
            raise ValueError("Precision needs to be either 'single' or 'double'.")

        self.max_iterations = max_iterations
        self.min_iterations = min_iterations

        # If V is not in a batch, put it in a batch of 1
        if len(V.shape) == 2:
            V = V[None, :, :]

        self._V = V.type(self._tensor_type)
        self._fix_neg = nn.Threshold(0., 1e-8)
        self._tolerance = tolerance
        self._prev_loss = None
        self._iter = 0
        self._test_conv = test_conv
        # self._gpu_id = gpu_id
        self._rank = rank
        self._generator = generator
        self._W, self._H = self._initialise_wh(init_method)

    def _initialise_wh(self, init_method):
        """
        Initialise basis and coefficient matrices according to `init_method`
        """
        if init_method == 'random':
            W = torch.unsqueeze(
                torch.from_numpy(self._generator.random((self._V.shape[1], self._rank), dtype=np.float64)), 0)
            H = torch.unsqueeze(
                torch.from_numpy(self._generator.random((self._rank, self._V.shape[2]), dtype=np.float64)), 0)
            if self._np_dtype is np.float32:
                W = W.float()
                H = H.float()
            return W, H

        elif init_method == 'nndsvd':
            W = np.zeros([self._V.shape[0], self._V.shape[1], self._rank])
            H = np.zeros([self._V.shape[0], self._rank, self._V.shape[2]])
            nv = nndsvd.Nndsvd()
            for i in range(self._V.shape[0]):
                vin = np.mat(self._V.cpu().numpy()[i])
                W[i, :, :], H[i, :, :] = nv.initialize(vin, self._rank, options={'flag': 0})

        elif init_method == 'nndsvda':
            W = np.zeros([self._V.shape[0], self._V.shape[1], self._rank])
            H = np.zeros([self._V.shape[0], self._rank, self._V.shape[2]])
            nv = nndsvd.Nndsvd()
            for i in range(self._V.shape[0]):
                vin = np.mat(self._V.cpu().numpy()[i])
                W[i, :, :], H[i, :, :] = nv.initialize(vin, self._rank, options={'flag': 1})

        elif init_method == 'nndsvdar':
            W = np.zeros([self._V.shape[0], self._V.shape[1], self._rank])
            H = np.zeros([self._V.shape[0], self._rank, self._V.shape[2]])
            nv = nndsvd.Nndsvd()
            for i in range(self._V.shape[0]):
                vin = np.mat(self._V.cpu().numpy()[i])
                W[i, :, :], H[i, :, :] = nv.initialize(vin, self._rank, options={'flag': 2})
        elif init_method == 'nndsvd_min':
            W = np.zeros([self._V.shape[0], self._V.shape[1], self._rank])
            H = np.zeros([self._V.shape[0], self._rank, self._V.shape[2]])
            nv = nndsvd.Nndsvd()
            for i in range(self._V.shape[0]):
                vin = np.mat(self._V.cpu().numpy()[i])
                w, h = nv.initialize(vin, self._rank, options={'flag': 2})
                min_X = np.min(vin[vin > 0])
                h[h <= min_X] = min_X
                w[w <= min_X] = min_X
                # W= np.expand_dims(W, axis=0)
                # H = np.expand_dims(H, axis=0)
                W[i, :, :] = w
                H[i, :, :] = h
        # W,H=initialize_nm(vin, nfactors, init=init, eps=1e-6,random_state=None)
        W = torch.from_numpy(W).type(self._tensor_type)
        H = torch.from_numpy(H).type(self._tensor_type)
        return W, H

    @property
    def reconstruction(self):
        return self.W @ self.H

    @property
    def W(self):
        return self._W

    @property
    def H(self):
        return self._H

    @property
    def conv(self):
        try:
            return self._conv
        except:
            return 0

    @property
    def _kl_loss(self):
        # calculate kl_loss in double precision for better convergence criteria
        return (self._V * (self._V / self.reconstruction).log()).sum(dtype=torch.float64) - self._V.sum(
            dtype=torch.float64) + self.reconstruction.sum(dtype=torch.float64)

    @property
    def generator(self):
        return self._generator

    @property
    def _loss_converged(self):
        """
        Check if loss has converged
        """
        if not self._iter:
            self._loss_init = self._kl_loss
        elif ((self._prev_loss - self._kl_loss) / self._loss_init) < self._tolerance:
            return True
        self._prev_loss = self._kl_loss
        return False

    def fit(self, beta=1):
        """
        Fit the basis (W) and coefficient (H) matrices to the input matrix (V) using multiplicative updates and
            beta divergence
        Args:
          beta: value to use for generalised beta divergence. Default is 1 for KL divergence
            beta == 2 => Euclidean updates
            beta == 1 => Generalised Kullback-Leibler updates
            beta == 0 => Itakura-Saito updates
        """
        with torch.no_grad():
            def stop_iterations():
                stop = (self._V.shape[0] == 1) and \
                       (self._iter % self._test_conv == 0) and \
                       self._loss_converged and \
                       (self._iter > self.min_iterations)
                if stop:
                    pass
                    # print("loss converged with {} iterations".format(self._iter))
                return [stop, self._iter]

            if beta == 2:
                for self._iter in range(self.max_iterations):
                    self._H = self.H * (self.W.transpose(1, 2) @ self._V) / (self.W.transpose(1, 2) @ (self.W @ self.H))
                    self._W = self.W * (self._V @ self.H.transpose(1, 2)) / (self.W @ (self.H @ self.H.transpose(1, 2)))
                    if stop_iterations()[0]:
                        self._conv = stop_iterations()[1]
                        break

            # Optimisations for the (common) beta=1 (KL) case.
            elif beta == 1:
                ones = torch.ones(self._V.shape).type(self._tensor_type)
                for self._iter in range(self.max_iterations):
                    print(self.W.shape, self._V.shape, self.H.shape)
                    ht = self.H.transpose(1, 2)
                    numerator = (self._V / (self.W @ self.H)) @ ht

                    denomenator = ones @ ht
                    self._W *= numerator / denomenator

                    wt = self.W.transpose(1, 2)
                    numerator = wt @ (self._V / (self.W @ self.H))
                    denomenator = wt @ ones
                    self._H *= numerator / denomenator
                    if stop_iterations()[0]:
                        self._conv = stop_iterations()[1]
                        break

            else:
                for self._iter in range(self.max_iterations):
                    self._H = self.H * ((self.W.transpose(1, 2) @ (((self.W @ self.H) ** (beta - 2)) * self._V)) /
                                        (self.W.transpose(1, 2) @ ((self.W @ self.H) ** (beta - 1))))
                    self._W = self.W * ((((self.W @ self.H) ** (beta - 2) * self._V) @ self.H.transpose(1, 2)) /
                                        (((self.W @ self.H) ** (beta - 1)) @ self.H.transpose(1, 2)))
                    if stop_iterations()[0]:
                        self._conv = stop_iterations()[1]
                        break


if __name__ == "__main__":
    X = torch.tensor(np.loadtxt(sys.argv[1], dtype=float, delimiter=",", skiprows=1)).T
    import numpy as np
    from sklearn.decomposition import NMF
    model = NMF(n_components=10, init='random', random_state=0, max_iter=1000)
    W = model.fit_transform(X)
    H = model.components_
    print(W.shape, H.shape)
    print(H.sum(axis=-1))
    # print(W)
#     print(data)
#     nmf = NMF(data, 3, max_iterations=5)
#     nmf.fit()
#     print(nmf.W)
