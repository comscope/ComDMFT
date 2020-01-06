import numpy as np
import h5py
from pyglib.estructure import dos as _dos
from future_builtins import zip


with h5py.File('GBANDS_0.h5', 'r') as f:
    k_start = f['/IKP_START'][0]
    k_end = f['/IKP_END'][0]
    e_kn = []
    psi_kn = []
    for k in range(k_start, k_end+1):
        e_kn.append(f['/ISPIN_1/IKP_{}/ek'.format(k)][...])
        psi_kn.append(f['/ISPIN_1/IKP_{}/ISYM_1/EVEC'.format(k)][...].T)

with h5py.File('BAREHAM_0.h5', 'r') as f:
    for k in range(k_start, k_end+1):
        psi_kn[k-1] = f['/IKP_{}/T_PSIK0_TO_HK0_BASIS'.format(k)][...].T. \
                dot(psi_kn[k-1])

with h5py.File('GLOG.h5', 'r') as f:
    efermi = f['/e_fermi'][0]

e_kn = np.asarray(e_kn) - efermi

with open('eigenvalues.dat', 'w') as f:
    for k in range(k_start, k_end+1):
        for n in range(e_kn.shape[1]):
            f.write(' {:10d} {:10d} {:20.14f}\n'.format(n+1,k,e_kn[k-1,n]))

with open('psi0_psig.dat', 'w') as f:
    for k in range(k_start, k_end+1):
        for n in range(psi_kn[k-1].shape[0]):
            f.write(' {:10d} {:10d}\n'.format(k, n+1))
            for m in range(psi_kn[k-1].shape[1]):
                f.write(' {:20.14f} {:20.14f}\n'.format( \
                        psi_kn[k-1][n, m].real, psi_kn[k-1][n, m].imag))

