import numpy as np
from builtins import range
import h5py

sz_list = []

with h5py.File("GPARAM.h5", 'r') as f:
    num_imp = f["/num_imp"][0]
    for i in range(num_imp):
        sz_list.append(f["/IMPURITY_"+str(i+1)+"/SZ"][...].T)

nks_list = []
for sz in sz_list:
    nks_list.append(np.zeros_like(sz))
    nks_list[-1][4,4] = nks_list[-1][6,6] = 1.0

with h5py.File("WH_RL_INIT.h5",'w') as f:
    for i,nks in enumerate(nks_list):
        f["/IMPURITY_"+str(i+1)+"/NKS"] = nks.T
