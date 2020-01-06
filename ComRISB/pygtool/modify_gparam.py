from __future__ import print_function
import numpy as np
import h5py


f = h5py.File('GPARAM.h5','a')

f['/gmaxiter'][...] = [1]
f['/giembeddiag'][...] = [-1]
f['/nval_top_ityp'][...] = [4]
f['/gimix'][...] = [-1]
f['/dc_mode'][...] = [1]

f.close()
