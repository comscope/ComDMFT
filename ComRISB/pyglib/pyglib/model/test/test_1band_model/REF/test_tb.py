import numpy as np
from ase.dft import kpoints
import pyglib.gutz.ginput as ginput
import pyglib.model.tbASE as tb

# The following is a simple test for the 1-d Hubbard model.

a = tb.AtomsTB("N", [(0, 0, 0)], cell=(1, 1, 1))
a.set_orbitals_spindeg()
aTB = tb.TB(a)
aTB.set_hop([
    ((1, 0, 0), 0, 0, -1),
    ((-1, 0, 0), 0, 0, -1),
    ((0, 0, 0), 0, 0,  0), ])

kps_size = (100, 1, 1)
kps = kpoints.monkhorst_pack(kps_size)

num_k = len(kps)
kps_wt = 1.0 / num_k * np.ones((num_k))
if aTB.Atoms.spindeg:
    kps_wt *= 2
num_e = 1.0
num_band_max = 1

# GPARAMBANDS.h5
h1e_list = [np.array([[0, 0], [0, 0]], dtype=np.complex)]
ginput.save_gparambands(kps_wt, num_e, num_band_max, h1e_list=h1e_list)

sigma_list = [np.identity(2, dtype=np.int32)]
v2e = np.zeros((2, 2, 2, 2), dtype=np.complex)
v2e[0, 0, 0, 0] = v2e[0, 0, 1, 1] = v2e[1, 1, 0, 0] = v2e[1, 1, 1, 1] = 6.0
sz = np.asarray(np.diag((1,-1)),dtype=np.complex)

# GPARAM.h5
ginput.save_gparam(sigma_list=sigma_list, iembeddiag=-1, v2e_list=[v2e],
        sz_list=[sz])

# BAREHAM_0.h5
aTB.save_bareham(kps)
