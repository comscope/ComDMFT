import h5py
import numpy as np

evals = np.loadtxt('eigenvalues.dat')
evals = evals.T
num = 22
numk = evals.shape[1]/num

with h5py.File('data.h5', 'w') as f:

    for ik, ibase in enumerate(range(0, evals.shape[1], num)):
        f['/evals/ik_{}'.format(ik)] = evals[2,ibase:ibase+num]

    for ik in range(numk):
        orbs = np.loadtxt('projector_{}.dat'.format(ik+1))
        orbs = orbs.T
        orbs = orbs[3] + 1.j*orbs[4]
        orbs = orbs.reshape((2, 22, 5))
        orbs = np.swapaxes(orbs, 1, 2)
        orbs = orbs.reshape((10,22))
        orbs = orbs.T.conj()   # <psi_k | loc orb>
        f['/psi_orb/ik_{}'.format(ik)] = orbs
