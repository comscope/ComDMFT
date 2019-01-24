import numpy as np
import h5py
import sys
from pyglib.math.matrix_basis import sigmatomatrixbasis


def save_gparam(iso=1, ispin=1, ityp_list=[0], imap_list=[0], na2_list=[2],
        imix=0, amix=0.2, nbreset=0, iembeddiag=0, sigma_list=None,
        v2e_list=None, max_iter=10000, sx_list=None, sy_list=None,
        sz_list=None, lx_list=None, ly_list=None, lz_list=None,
        utrans_list=None, ldc=0, u_avg_list=None, j_avg_list=None,
        nelf_list=None, rotations_list=None, lie_odd_params_list=None,
        lie_even_params_list=None, jgenerator_list=None,
        sp_rotations_list=None,
        nval_bot_list=[0], nval_top_list=[2]):

    if '-n' in sys.argv:
        max_iter = np.int32(sys.argv[sys.argv.index('-n') + 1])

    with h5py.File('GPARAM.h5', 'w') as f:
        f['/iso'] = np.asarray([iso], dtype=np.int32)
        f['/ispin'] = np.asarray([ispin], dtype=np.int32)
        f['/num_imp'] = np.asarray([len(ityp_list)], dtype=np.int32)
        f['/ITYP_IMP'] = np.asarray(ityp_list, dtype=np.int32) + 1
        f['/IMAP_IMP'] = np.asarray(imap_list, dtype=np.int32) + 1
        f['/na2_imp'] = np.asarray(na2_list, dtype=np.int32)
        f['/gimix'] = np.asarray([imix], dtype=np.int32)
        f['/gamix'] = [amix]
        f['/gnbreset'] = np.asarray([nbreset], dtype=np.int32)
        f['/giembeddiag'] = np.asarray([iembeddiag], dtype=np.int32)
        f['/gmaxiter'] = np.asarray([max_iter], dtype=np.int32)
        f['/dc_mode'] = np.asarray([ldc], dtype=np.int32)

        f['/dc_u_avg'] = u_avg_list
        f['/dc_j_avg'] = j_avg_list
        if nelf_list is not None:
            f['/dc_nelf_list'] = nelf_list

        if nval_bot_list is not None:
            f['/nval_bot_ityp'] = np.asarray(nval_bot_list)
            f['/nval_top_ityp'] = np.asarray(nval_top_list)

        # for sigma_list, v2e_list
        dim_hs_list = []
        for i, na2 in enumerate(na2_list):
            if sigma_list is None:
                sigma = np.arange(1,na2**2+1,dtype=np.int32).reshape((na2,na2))
            else:
                sigma = sigma_list[i]

            f['/IMPURITY_'+str(i+1)+'/SIGMA_STRUCT'] = sigma.T
            hermitian_mb = sigmatomatrixbasis(sigma)
            dim_hs_list.append(len(hermitian_mb))
            f['/IMPURITY_'+str(i+1)+'/HS'] = np.swapaxes(hermitian_mb,1,2)

            if v2e_list is None:
                v2e = np.zeros((na2,na2,na2,na2), dtype=np.complex)
            else:
                v2e = v2e_list[i]
            f['/IMPURITY_'+str(i+1)+'/V2E'] = v2e.T

            if sx_list is not None:
                f['/IMPURITY_'+str(i+1)+'/SX'] = sx_list[i].T
            if sy_list is not None:
                f['/IMPURITY_'+str(i+1)+'/SY'] = sy_list[i].T
            if sz_list is not None:
                f['/IMPURITY_'+str(i+1)+'/SZ'] = sz_list[i].T
            if lx_list is not None:
                f['/IMPURITY_'+str(i+1)+'/LX'] = lx_list[i].T
            if ly_list is not None:
                f['/IMPURITY_'+str(i+1)+'/LY'] = ly_list[i].T
            if lz_list is not None:
                f['/IMPURITY_'+str(i+1)+'/LZ'] = lz_list[i].T
            if utrans_list is not None:
                f['/IMPURITY_'+str(i+1)+'/DB_TO_SAB']\
                        = utrans_list[i].T
            if rotations_list is not None:
                f['/IMPURITY_'+str(i+1)+'/rotations'] = rotations_list[i]
            if lie_odd_params_list is not None:
                f['/IMPURITY_'+str(i+1)+'/lie_odd_params'] = \
                        lie_odd_params_list[i]
                f['/IMPURITY_'+str(i+1)+'/lie_even_params'] = \
                        lie_even_params_list[i]
                f['/IMPURITY_'+str(i+1)+'/JGENERATOR'] = \
                        np.swapaxes(jgenerator_list[i],1,2)
                f['/IMPURITY_'+str(i+1)+'/SP_ROTATIONS'] = \
                        np.swapaxes(sp_rotations_list[i],1,2)
                f['/IMPURITY_'+str(i+1)+'/nsym_odd'] = \
                        np.asarray(sp_rotations_list[i]).shape[0:1]

        f['/dim_hs_imp'] = np.asarray(dim_hs_list, dtype=np.int32)


def save_gparambands(kpt_wt, nelectron, nbmax, iso=1, ispin=1,
        symop=1, symie=1, ensemble=0,
        ismear=0, delta=0.01, ne_list=None, h1e_list=None):
    with h5py.File('GPARAMBANDS.h5', 'w') as f:
        f['/iso'] = np.asarray([iso], dtype=np.int32)
        f['/ispin'] = np.asarray([ispin], dtype=np.int32)
        f['/kptdim'] = np.asarray(list(kpt_wt.shape), dtype=np.int32)
        f['/nbmax'] = np.asarray([nbmax], dtype=np.int32)
        if ne_list is not None:
            assert ne_list.shape[0] == kpt_wt.shape[0], \
                ' error: ne_list and kpt_wt do not match!'
            f['/NE_LIST'] = np.asarray(ne_list, dtype=np.int32)  # shape (nk,3)
        f['/kptwt'] = kpt_wt
        f['/ismear'] = np.asarray([ismear], dtype=np.int32)
        f['/ensemble'] = [ensemble]
        f['/delta'] = [delta]
        f['/nelectron'] = [nelectron]
        f['/symnop'] = np.asarray([symop], dtype=np.int32)
        f['/symie'] = np.asarray([symie], dtype=np.int32)
        if h1e_list is None:
            h1e_list = [np.zeros((2, 2), dtype=np.complex)]
        for i, h1e in enumerate(h1e_list):
            # Fortran convention.
            f['/IMPURITY_' + str(i + 1) + '/H1E'] = h1e.T
