import numpy as np
from ase.dft import kpoints
import pyglib.gutz.ginput as ginput
import pyglib.model.tbASE as tb
import scipy.constants as s_const

# The following is a simple test for the 1-d Hubbard model.

a = tb.AtomsTB("N", [(0, 0, 0)], cell=(1, 1, 1))
a.set_orbitals_spindeg(orbitals=[("s","p")])
aTB = tb.TB(a)
aTB.set_hop([
        (( 1,0,0),0,0,-1),
        ((-1,0,0),0,0,-1),
        (( 1,0,0),1,1,-1),
        ((-1,0,0),1,1,-1),
        ((0,0,0),0,1,0.25),
        ((0,0,0),1,0,0.25)])

kps_size = (100, 1, 1)
kps = kpoints.monkhorst_pack(kps_size)

num_k = len(kps)
kps_wt = 1.0 / num_k * np.ones((num_k))
if aTB.Atoms.spindeg:
    kps_wt *= 2
num_e = 1.2
num_band_max = 2

# GPARAMBANDS.h5
h1e_list = [np.array([[0,0,0.25,0],[0,0,0,0.25],[0.25,0,0,0],[0,0.25,0,0]], \
        dtype=np.complex)]
ginput.save_gparambands(kps_wt, num_e, num_band_max, h1e_list=h1e_list)
# GPARAM.h5
sigma = np.zeros((4,4), dtype=np.int32)
sigma[0::2, 0::2] = np.arange(1,5).reshape((2,2))
sigma[1::2, 1::2] = sigma[0::2, 0::2]
v2e = np.zeros((4, 4, 4, 4), dtype=np.complex)
v2e[0, 0, 0, 0] = v2e[0, 0, 1, 1] = v2e[1, 1, 0, 0] = v2e[1, 1, 1, 1] = 8.0
v2e[2:, 2:, 2:, 2:] = v2e[0:2, 0:2, 0:2, 0:2]
ryd2ev = s_const.physical_constants['Rydberg constant times hc in eV'][0]
sz = np.asarray(np.diag((1,-1,1,-1)),dtype=np.complex)

ginput.save_gparam(na2_list=[4], iembeddiag=-1, imix=-1, sigma_list=[sigma],
        v2e_list=[v2e], sz_list=[sz], nval_top_list=[4])

# BAREHAM_0.h5
aTB.save_bareham(kps)
