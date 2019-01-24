from __future__ import print_function
import h5py


'''Directly modify the Coulomb U-matrix of impurity assuming Slater-Condon type
parametrization.

Files needed:
    GPARAM.h5, EMBED_HAMIL_imp.h5
'''


# impurity index
imp = 1

# U/J parameters to be specified (eV).
u = 3.0
j = 0.7

from pyglib.mbody.h5coulomb_matrix import get_v2e

v2e = get_v2e(u, j, imp)
with h5py.File('EMBED_HAMIL_{}.h5'.format(imp)) as f:
    f['/V2E'][...] = v2e.T
