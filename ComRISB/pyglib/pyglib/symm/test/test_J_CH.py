from __future__ import print_function

import numpy as np
import unittest
import pyglib.symm.angular_momentum_1p as am
from builtins import zip

J = am.get_J_vector([3], 'RH')
S = am.get_S_vector_JJ([3])
L = am.get_L_vector_JJ([3])

class KnowValues(unittest.TestCase):
    def test_BasisTransformation(self):
        for _J, _S, _L in zip(J, S, L):
            self.assertTrue(np.max(np.abs(_J-_S-_L))<1.e-6)


if __name__ == "__main__":
    print("Tests for J_CH_Transformation.")
    unittest.main()

