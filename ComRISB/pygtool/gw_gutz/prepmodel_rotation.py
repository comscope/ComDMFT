import numpy as np
import h5py


def get_lowdin_orthonorm(a):
    assert a.shape[1] <= a.shape[0], 'Error: n_orb > n_psi!'

    res = np.zeros((a.shape[0],a.shape[0]), dtype=a.dtype)
    u, _, vh = np.linalg.svd(a)
    res[:, 0:vh.shape[0]] = u[:, 0:vh.shape[0]].dot(vh)

    if vh.shape[0] < res.shape[0]:
        b = res[:, 0:vh.shape[0]].dot(res[:, 0:vh.shape[0]].T.conj())
        _, v = np.linalg.eigh(b)
        res[:, vh.shape[0]:] = v[:, 0:res.shape[0]-vh.shape[0]]

    return res


def get_hamilt_from_ev(w, v):
    vh = v.T.conj()
    for i in range(vh.shape[0]):
        vh[i] *= w
    return vh.dot(v)


def h5set_model():
    numk = 396
    wk = 1./numk
    numorb = 10
    nbmax = 22
    h1e = np.zeros((numorb, numorb), dtype=complex)

    # Unitary transformation from {xy, yz, z2, zx, x2-y2} to complex Harmonics.
    ucm2cu = np.zeros((5, 5), dtype=complex)
    ucm2cu[0, 0] =  1.j/np.sqrt(2.); ucm2cu[4, 0] = -1.j/np.sqrt(2.)
    ucm2cu[1, 1] = -1.j/np.sqrt(2.); ucm2cu[3, 1] = -1.j/np.sqrt(2.)
    ucm2cu[2, 2] = 1.0
    ucm2cu[1, 3] = 1./np.sqrt(2); ucm2cu[3, 3] = -1./np.sqrt(2.)
    ucm2cu[0, 4] = 1./np.sqrt(2.); ucm2cu[4, 4] = 1./np.sqrt(2.)

    from scipy.linalg import block_diag
    ucm2cu = block_diag(ucm2cu, ucm2cu)

    fm = h5py.File('BAREHAM_0.h5', 'w')
    with h5py.File('data.h5', 'r') as f:
        for ik in range(numk):
            evals = f['/evals/ik_{}'.format(ik)][...]
            fm['/IKP_{}/ek0'.format(ik+1)] = evals

            psi_orb = f['/psi_orb/ik_{}'.format(ik)][...]
            evecs = get_lowdin_orthonorm(psi_orb)

            # fix phase convension.
            evecs[:,numorb/2+1] *= -1.
            evecs[:,numorb/2+2] *= -1.
            evecs[:,numorb/2+4] *= -1.

            # Additional rotation:
            # Rotate from {xy, yz, z2, zx, x2-y2} to complex Harmonics.
            evecs[:, 0:numorb] = evecs[:, 0:numorb].dot(ucm2cu.T.conj())

            fm['/IKP_{}/T_PSIK0_TO_HK0_BASIS'.format(ik+1)] = \
                    evecs.T

            hmat = get_hamilt_from_ev(evals, evecs)
            fm['/IKP_{}/ISYM_1/HK0'.format(ik+1)] = hmat.T
            h1e += hmat[0:numorb, 0:numorb]*wk
    fm.close()

    assert np.allclose(h1e[0:numorb/2, 0:numorb/2], \
            h1e[numorb/2:numorb, numorb/2:numorb]), 'Check h1e!'

    with h5py.File('GPARAMBANDS.h5', 'w') as f:
        f['/iso/'] = [1]
        f['/ispin'] = [1]
        f['/kptdim'] = [numk]
        f['/nbmax'] = [nbmax]
        f['/NE_LIST'] = [[nbmax, 1, nbmax] for k in range(numk)]
        f['/kptwt'] = [wk for k in range(numk)]
        f['/ismear'] = [0]
        f['/delta'] = [0.01]
        f['/nelectron'] = [28.0]
        f['/symnop'] = [1]
        f['/symie'] = [1]
        f['/IMPURITY_1/H1E'] = block_diag(h1e[0:numorb/2, 0:numorb/2],
                h1e[0:numorb/2, 0:numorb/2]).T
        f['/IMPURITY_2/H1E'] = block_diag(h1e[numorb/2:numorb, numorb/2:numorb],
                h1e[numorb/2:numorb, numorb/2:numorb]).T


if __name__ == '__main__':
    h5set_model()
