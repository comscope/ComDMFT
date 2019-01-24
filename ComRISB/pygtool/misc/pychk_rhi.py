import h5py, numpy

f1 = h5py.File('EMBED_HAMIL_PHIMAT_1.h5', 'r')
f2 = h5py.File('EMBED_HAMIL_ANALYSIS_1.h5', 'r')
for i in range(11):
    path = '/valence_block_{}/PHI'.format(i)
    if path not in f1:
        continue
    phi = f1['/valence_block_{}/PHI'.format(i)][()].T
    rho = f2['/valence_block_{}/RHO'.format(i)][()].T
    print numpy.max(numpy.abs(rho - phi.dot(phi.T.conj()))), 'right'
    print numpy.max(numpy.abs(rho - phi.T.conj().dot(phi))), 'wrong'
