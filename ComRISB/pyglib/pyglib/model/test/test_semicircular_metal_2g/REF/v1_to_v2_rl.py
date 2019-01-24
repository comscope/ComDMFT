import h5py
import numpy as np


np.set_printoptions(precision=4,suppress=True, linewidth=100)

with h5py.File('glog.h5', 'r') as f:
    rmat = f['/Impurity_1/GA_R'][:3,:3]
    lmat = f['/Impurity_1/GA_La'][:3,:3]

umat = np.array([
    [0, 1, 0],
    [1, 0, 0],
    [0, 0, 1]], dtype=complex)

rmat = umat.T.conj().dot(rmat.dot(umat))
lmat = umat.T.conj().dot(lmat.dot(umat))

rmat2 = np.zeros((6, 6), dtype=complex)
rmat2[0::2, 0::2] = rmat
rmat2[1::2, 1::2] = rmat2[0::2, 0::2]

lmat2 = np.zeros((6, 6), dtype=complex)
lmat2[0::2, 0::2] = lmat
lmat2[1::2, 1::2] = lmat2[0::2, 0::2]

with h5py.File('WH_RL_INIT.h5', 'w') as f:
    f['/IMPURITY_1/LAMBDA'] = lmat2
    f['/IMPURITY_1/R'] = rmat2
