from __future__ import print_function

import os
import unittest
import numpy as np
from scipy.linalg import expm
import h5py
from pyglib.math.matrix_util import yield_derivative_f_matrix


class KnowValues(unittest.TestCase):
    def test_derivative_exp_ix(self):
        # random Hermitian matrix of dimention 5.
        n = 5
        x = np.random.rand(n, n) + np.random.rand(n, n)*1.j
        x += np.conj(x.T)
        # Numerically calculate the derivative w.r.t. h.
        h = np.array([ \
            [0,              0, 1./np.sqrt(2.), 0, 0],
            [0,              0, 0,              0, 0],
            [1./np.sqrt(2.), 0, 0,              0, 0],
            [0,              0, 0,              0, 0],
            [0,              0, 0,              0, 0]])

        expix0 = expm(-1.j*x)
        delta = 1.e-7
        expix1 = expm(-1.j*(x + delta*h))
        partial0 = (expix1 - expix0)/delta
        f = lambda x: np.exp(-1.j*x)
        fp = lambda x: -1.j*np.exp(-1.j*x)
        partials = yield_derivative_f_matrix(x, [h], f, fp)
        for partial in partials:
            err = np.abs(np.max(partial - partial0))
            self.assertAlmostEqual(err, 0.)

    def test_derivative_tr_exp_ix_a(self):
        # random Hermitian matrix of dimention 5.
        n = 5
        a = np.random.rand(n, n)
        x = np.random.rand(n, n) + np.random.rand(n, n)*1.j
        x += np.conj(x.T)
        # Numerically calculate the derivative w.r.t. h.
        h = np.array([ \
            [0,              0, 1./np.sqrt(2.), 0, 0],
            [0,              0, 0,              0, 0],
            [1./np.sqrt(2.), 0, 0,              0, 0],
            [0,              0, 0,              0, 0],
            [0,              0, 0,              0, 0]])

        expix0 = expm(-1.j*x)
        delta = 1.e-7
        expix1 = expm(-1.j*(x + delta*h))
        partial0 = (np.trace(np.dot(expix1, a)) -
                np.trace(np.dot(expix0, a)))/delta
        f = lambda x: np.exp(-1.j*x)
        fp = lambda x: -1.j*np.exp(-1.j*x)
        partials = yield_derivative_f_matrix(x, [h], f, fp)
        for partial in partials:
            # trace
            _partial = np.sum(partial.T*a)
            err = np.abs(_partial - partial0)
            self.assertAlmostEqual(err, 0.)


if __name__=="__main__":
    print("Tests for matrix_util.")
    unittest.main()
