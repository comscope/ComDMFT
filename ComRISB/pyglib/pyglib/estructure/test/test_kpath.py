from __future__ import print_function
import numpy
from pyglib.estructure.gwannier import get_symkpath


kpath = get_symkpath()
print("lattice")
print(kpath._structure.lattice.matrix)
print(kpath._prim.lattice.matrix)

print("reciprocal lattice")
print(kpath._structure.lattice.reciprocal_lattice.matrix)
print(kpath._prim.lattice.reciprocal_lattice.matrix)

print("k_path")
print(kpath.kpath)
