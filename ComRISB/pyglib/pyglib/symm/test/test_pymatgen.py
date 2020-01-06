from __future__ import print_function

import unittest
import pickle

with open('mol.pickle', 'rb') as f:
    mol = pickle.load(f)

from pymatgen.symmetry.analyzer import PointGroupAnalyzer
analyzer = PointGroupAnalyzer(mol)
print(" sch_symbol = ", analyzer.sch_symbol)

from pymatgen.symmetry.analyzer import generate_full_symmops
symmops = generate_full_symmops(analyzer.symmops, tol=1.e-7)
print(" num symmops = ", len(symmops))


class KnowValues(unittest.TestCase):
    def test_PointGroupAnalyzer(self):
        self.assertTrue(analyzer.sch_symbol == 'D4h')
        self.assertTrue(len(symmops) == 16)


if __name__ == "__main__":
    print("Tests for pymatgen.")
    unittest.main()

