import h5py
import numpy as np

case = 'UO2_6'

dmat = np.zeros((14,14),dtype=complex)
with open('UO2_6.dmatup', 'r') as f:
    data = []
    for line in f.readlines()[2:]:
        data += [float(x) for x in line.split()]
    data = np.asarray(data)
    data = data[::2] + data[1::2]*1.j
    data = data.reshape((7,7))
    dmat[:7,:7] = data

with open('UO2_6.dmatdn', 'r') as f:
    data = []
    for line in f.readlines()[2:]:
        data += [float(x) for x in line.split()]
    data = np.asarray(data)
    data = data[::2] + data[1::2]*1.j
    data = data.reshape((7,7))
    dmat[7:,7:] = data

with open('UO2_6.dmatud', 'r') as f:
    data = []
    for line in f.readlines()[2:]:
        data += [float(x) for x in line.split()]
    data = np.asarray(data)
    data = data[::2] + data[1::2]*1.j
    data = data.reshape((7,7))
    dmat[:7,7:] = data
    dmat[7:,:7] = data.T.conj()

print(np.trace(dmat))

with h5py.File('dmat.h5', 'w') as f:
    f['/dmat'] = dmat
