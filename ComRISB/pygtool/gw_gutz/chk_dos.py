import numpy as np
import h5py
from pyglib.estructure import dos as _dos
from future_builtins import zip


def count_elecn(dos):
    ne = 0.
    for e_kn in dos.e_skn:
        for w, e_n in zip(dos.w_k, e_kn):
            for e in e_n:
                if e > 0.0:
                    continue
                ne += w
    return ne


with h5py.File("GPARAMBANDS.h5", 'r') as f:
    w_k = f["/kptwt"][...]
    bnd_ne = f["/NE_LIST"][...]
    nsym = f["/symnop"][0]
    nkp = f["/kptdim"][0]
    nbmax = f["/nbmax"][0]

with h5py.File('BAREHAM_0.h5', 'r') as f:
    e_kn = []
    for k in range(nkp):
        e_kn.append(f['/IKP_{}/ek0'.format(k+1)][...])
e_skn = np.asarray([e_kn])

with h5py.File('data.h5', 'r') as f:
    psi_knf = []
    for k in range(nkp):
        tmp = f['/psi_orb/ik_{}'.format(k)][...]
        psi_knf.append([tmp[n,:].dot(tmp[n,:].conj()) for n in range(nbmax)])
psi_sknf = [psi_knf]

window = (e_skn.min(), e_skn.max())
dos = _dos.DOS(w_k, e_skn,  width=0.05, window=window, npts=1001)
energies = dos.get_energies()
dos_t = dos.get_dos_t()
dos_f = dos.get_dos_component(psi_sknf)

print 'num elecn = {}'.format(count_elecn(dos))

import matplotlib.pyplot as plt
plt.figure(0)
plt.plot(energies, dos_t[0].real, '-', label='t')
plt.plot(energies, dos_f[0].real, 'o', label='f')
plt.legend()
plt.xlim(-0.5,0.5)
plt.xlabel("E (eV)")
plt.ylabel("DOS (states/f.u.)")
plt.show()
