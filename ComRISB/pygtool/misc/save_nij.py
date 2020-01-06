import h5py
from pyglib.io.h5io import get_coo_matrix, write_coo_matrix

fsrc = h5py.File('EMBED_HAMIL_ANALYSIS_1.h5', 'r')
fdst = h5py.File('mbody.h5', 'a')

for val in range(1,10):
    for i in range(10):
        for j in range(i+1):
            n = get_coo_matrix(fsrc, '/valence_block_{}/NP_{}_{}'.format(\
                    val, i+1, j+1))
            path = '/d/valence_block_{}/n_{}_{}'.format(val,i,j)
            if path in fdst:
                del fdst[path]
            write_coo_matrix(fdst, path, n)
