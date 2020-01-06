#!/usr/bin/env python
from __future__ import print_function
import h5py, numpy, os
from pyglib.math.matrix_util import sym_array
from pyglib.gutz.usrqa import get_usr_input
from pyglib.basic.units import Ryd_eV
import pyglib.mbody.double_counting as dc


def get_b_field_list():
    '''get the list of magnetic field by q&a.
    '''
    with h5py.File('ginit.h5', 'r') as f:
        unit = f['/usrqa/unit'][()]
        unique_magmom_direction_list = \
                f["/usrqa/unique_magmom_direction_list"][()]
    if 'ryd' in unit:
        unit = Ryd_eV
    else:
        unit = 1.0

    with h5py.File('GPARAM.h5', 'r') as f:
        ispin = f["/ispin"][0]
        if ispin == 1:
            return None, None
        num_imp = f["/num_imp"][0]
        imap_imp = f['/IMAP_IMP'][()]

    print(' total {} impurities with equivalence indices \n {}'.format( \
            num_imp, imap_imp))
    bvec_list = []
    imp_unique = 0
    for imp in range(num_imp):
        if imap_imp[imp] == imp + 1:
            print('\n IMPURITY {}'.format(imp+1))
            vec = unique_magmom_direction_list[imp_unique]
            print("\n magnetic moment direction is {}".format(vec) +
                    "\n in global coordinate system.")
            bm = raw_input(
                    '\n please enter the magnitude of local B-field:'+\
                    '\n in unit of eV per Bohr magneton: \n ')
            bm = float(bm)/unit
            bvec_list.append(bm*vec)
        else:
            bvec_list.append(bvec_list[imap_imp[imp] - 1])
    givext = get_usr_input( \
            '\n Is the external field applied only at initial step (0) \n'+ \
            ' or fixed through all the iterations (1).', ['0', '1'])
    givext = int(givext)
    return bvec_list, givext


def get_sym_2darray(a, imp=1):
    '''get symmetrized 2d array with error reporting.
    '''
    with h5py.File('GPARAM.h5', 'r') as f:
        path = '/IMPURITY_{}/SP_ROTATIONS'.format(imp)
        if path in f:
            r_list = numpy.swapaxes(f[path][()], 1, 2)
            a_sym = sym_array(a, r_list)
        else:
            a_sym = a
        sym_err = numpy.max(numpy.abs(a-a_sym))
    return a_sym, sym_err


def get_vext_list(bvec_list):
    '''get the list of magnetic potential matrix in the single particle space,
    given the list of local magnetic field.
    '''
    vext_list = []
    max_sym_err = 0.
    with h5py.File('GPARAM.h5', 'r') as f:
        for imp, bvec in enumerate(bvec_list):
            prepath = "/IMPURITY_" + str(imp + 1)
            sx = f[prepath+'/SX'][()].T
            sy = f[prepath+'/SY'][()].T
            sz = f[prepath+'/SZ'][()].T
            vext = (bvec[0]*sx + bvec[1]*sy + bvec[2]*sz)*2
            lx = f[prepath+'/LX'][()].T
            ly = f[prepath+'/LY'][()].T
            lz = f[prepath+'/LZ'][()].T
            vext += bvec[0]*lx + bvec[1]*ly + bvec[2]*lz
            vext_sym, sym_err = get_sym_2darray(vext, imp+1)
            max_sym_err = max(max_sym_err, sym_err)
            vext_list.append(vext_sym)
    if max_sym_err > 1.e-5:
        print(" Warning:")
    print(' maximal symmetrization error of vext = {}'.format(max_sym_err))
    return vext_list


def get_vext_given_1pdm_list(dm_list):
    '''get the external potential for lda+u calculation
    given the initial one-particle density matrix.
    '''
    with h5py.File("GPARAM.h5", "r") as f:
        javg_list = f["/dc_j_avg"][()]
        uavg_list = f["/dc_u_avg"][()]
        v2e_list = []
        for i in range(f["/num_imp"][0]):
            v2e_list.append(f["/IMPURITY_{}/V2E".format(i+1)][()].T)

    vext_list = [-dc.get_vdc_hf(v2e, dm)+ \
            dc.get_vdc_fll(uavg, javg, numpy.trace(dm))*\
            numpy.eye(dm.shape[0]) \
            for v2e, dm, uavg, javg in zip(v2e_list, dm_list, \
            uavg_list, javg_list)]
    return vext_list


def chk_local_one_body(vext_list):
    if not os.path.isfile("GPARAMBANDS.h5"):
        return
    db2sab_list = []
    sx_list = []
    sy_list = []
    sz_list = []
    lx_list = []
    ly_list = []
    lz_list = []
    with h5py.File('GPARAM.h5', 'r') as f:
        num_imp = f["/num_imp"][0]
        for imp in range(num_imp):
            prepath = "/IMPURITY_" + str(imp + 1)
            db2sab_list.append(f[prepath+"/DB_TO_SAB".format(imp+1)] \
                    [()].T)
            sx_list.append(f[prepath+"/SX"][()].T)
            sy_list.append(f[prepath+"/SY"][()].T)
            sz_list.append(f[prepath+"/SZ"][()].T)
            lx_list.append(f[prepath+"/LX"][()].T)
            ly_list.append(f[prepath+"/LY"][()].T)
            lz_list.append(f[prepath+"/LZ"][()].T)
    h1e_list = []
    with h5py.File("GPARAMBANDS.h5", 'r') as f:
        for imp in range(num_imp):
            h1e = f["/IMPURITY_{}/H1E".format(imp+1)][()].T
            db2sab = db2sab_list[imp]
            h1e = db2sab.T.conj().dot(h1e).dot(db2sab)
            h1e_list.append(h1e)
    # start actually analysis
    numpy.set_printoptions(precision=2, suppress=True)
    for imp in range(num_imp):
        h = h1e_list[imp] + vext_list[imp]
        w, v = numpy.linalg.eigh(h)
        dm = v[:,0:1].dot(v[:,0:1].T.conj())
        dm_sym, _ = get_sym_2darray(dm, imp+1)
        print(" <Sx> = {}".format(numpy.sum(dm_sym*sx_list[imp])))
        print(" <Sy> = {}".format(numpy.sum(dm_sym*sy_list[imp])))
        print(" <Sz> = {}".format(numpy.sum(dm_sym*sz_list[imp])))
        print(" <Lx> = {}".format(numpy.sum(dm_sym*lx_list[imp])))
        print(" <Ly> = {}".format(numpy.sum(dm_sym*ly_list[imp])))
        print(" <Lz> = {}".format(numpy.sum(dm_sym*lz_list[imp])))


def h5wrt_gmagnet(vext_list, g_ivext, fn='GVEXT.h5'):
    with h5py.File(fn, 'w') as f:
        for imp, vext in enumerate(vext_list):
            f['/IMPURITY_{}/VEXT'.format(imp+1)] = vext.T
        f['/givext'] = [g_ivext]


def init_magnet_conf():
    '''
    initialize the the magnetic configuration for magnetic calculation.
    '''
    bvec_list, g_ivext = get_b_field_list()
    if bvec_list is None:
        return
    vext_list = get_vext_list(bvec_list)
    chk_local_one_body(vext_list)
    h5wrt_gmagnet(vext_list, g_ivext)


def init_magnet_conf_with_init_dm(dm_list):
    '''
    initialize the the magnetic configuration for magnetic calculation
    based on given initial one-particle density-matrix.
    '''
    vext_list = get_vext_given_1pdm_list(dm_list)
    chk_local_one_body(vext_list)
    h5wrt_gmagnet(vext_list, g_ivext=0)


if __name__ == "__main__":
    init_magnet_conf()
