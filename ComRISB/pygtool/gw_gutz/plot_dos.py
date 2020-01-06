import numpy as np
import h5py
from pyglib.estructure import dos as _dos
from future_builtins import zip


with h5py.File("../GPARAMBANDS.h5", 'r') as f:
    w_k = f["/kptwt"][...]
    bnd_ne = f["/NE_LIST"][...]
    nsym = f["/symnop"][0]
    nkp = f["/kptdim"][0]
    nbmax = f["/nbmax"][0]

with h5py.File('GLOG.h5', 'r') as f:
    e_fermi = f['/e_fermi'][0]

with h5py.File('GBANDS_0.h5', 'r') as f:
    k_start = f['/IKP_START'][0]
    k_end = f['/IKP_END'][0]
    e_kn = []
    psi_knf = []
    for k in range(k_start, k_end+1):
        e_kn.append(f['/ISPIN_1/IKP_{}/ek'.format(k)][...]-e_fermi)
        tmp = f['/ISPIN_1/IKP_{}/ISYM_1/EVEC'.format(k)][:,0:10]
        psi_knf.append([tmp[n,:].dot(tmp[n,:].conj()) \
                for n in range(nbmax)])
e_skn = np.asarray([e_kn])
psi_sknf = [psi_knf]

window = (e_skn.min(), e_skn.max())
dos = _dos.DOS(w_k, e_skn,  width=0.05, window=window, npts=1001)
energies = dos.get_energies()
dos_t = dos.get_dos_t()
dos_f = dos.get_dos_component(psi_sknf)

import matplotlib.pyplot as plt
plt.figure(0)
plt.plot(energies, dos_t[0].real, '-', label='t')
plt.plot(energies, dos_f[0].real, 'o', label='f')
plt.legend()
plt.xlim(-0.5,0.5)
plt.xlabel("E (eV)")
plt.ylabel("DOS (states/f.u.)")
plt.show()
