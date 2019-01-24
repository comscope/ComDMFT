'''Directly modify the Coulomb U-matrix of impurity
assuming Slater-Condon type parametrization.

Files needed:
    GPARAM.h5.
'''
import h5py


def get_v2e(u, j, imp):
    ryd2ev = 13.605698065894
    with h5py.File('GPARAM.h5', 'r') as f:
        iso = f['/iso'][0]
        # Upper-case Fortran convention to c-convention
        utrans = f['/IMPURITY_{}/DB_TO_SAB'.format(imp)][...].T

    # Get l quantum number
    l = (utrans.shape[0]/2-1)/2

    if iso == 2:
        from pyglib.math.matrix_util import trans_JJ_to_CH_sup_sdn
        u_cmplx_harm_to_rel_harm = trans_JJ_to_CH_sup_sdn(l).T
        utrans = u_cmplx_harm_to_rel_harm.dot(utrans)

    from pyglib.mbody.coulomb_matrix import U_matrix
    v2e, u_avg, j_avg = U_matrix('slater-condon', l, U_int=u,
            J_hund=j, T=utrans)
    print(' u_avg = {} j_avg = {}'.format(u_avg, j_avg))

    # CyGutz works with Rydberg units.
    v2e /= ryd2ev
    return v2e

if __name__ =='__main__':
    # impurity index
    imp = 1

    # U/J parameters to be specified (eV).
    u = 6.0
    j = 0.7

    v2e = get_v2e(u, j, imp)
    with h5py.File('EMBED_HAMIL_{}.h5'.format(imp)) as f:
        f['/V2E'][...] = v2e.T
