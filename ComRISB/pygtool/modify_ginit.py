from __future__ import print_function
import numpy as np
import h5py


with h5py.File('ginit.h5','a') as f:

    f['/usrqa/unique_j_list_ev'][()] *= 1.0
    f['/usrqa/unique_u_list_ev'][()] *= 6.0
    f['/usrqa/ldc'][()] = 1

