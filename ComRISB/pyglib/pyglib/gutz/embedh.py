from __future__ import print_function
# analyze the embeding hamiltonian
import os, h5py, numpy, shutil
from scipy.linalg import block_diag
from pyglib.math.matrix_util import \
        get_utrans_orbfast_supdn_to_spin_fast_supdn
from pyglib.math.matrix_util import unitary_transform_coulomb_matrix


def twrite_2d_array_real(a, fname, shreshold=1.e-10):
    '''write 2d array (real) in text form of i, j, a_ij.
    '''
    with open(fname, 'w') as f:
        for i, a1 in enumerate(a):
            for j, a12 in enumerate(a1):
                if numpy.abs(a12) > shreshold:
                    f.write('{:2d} {:2d} {:16.12f}\n'.format(\
                            i+1, j+1, a12.real))


def twrite_2d_array_cmplx(a, fname, shreshold=1.e-10):
    '''write 2d array (real) in text form of i, j, a_ij.
    '''
    with open(fname, 'w') as f:
        for i, a1 in enumerate(a):
            for j, a12 in enumerate(a1):
                if numpy.abs(a12) > shreshold:
                    f.write('{:2d} {:2d} {:16.12f} {:16.12f}\n'.format(\
                            i+1, j+1, a12.real, a12.imag))


def twrite_u_matrix_real(v2e, fname, shreshold=1.e-10):
    '''write u-matrix in text form of i, k, l, j, v2e_(ij)(kl).
    '''
    with open(fname, 'w') as f:
        for i, v1 in enumerate(v2e):
            for j, v12 in enumerate(v1):
                for k, v123 in enumerate(v12):
                    if i == k:
                        continue
                    for l, v1234 in enumerate(v123):
                        if j == l:
                            continue
                        if numpy.abs(v1234) > shreshold:
                            # order of cd_i cd_k c_l c_j
                            # absorbing 1/2 factor
                            f.write('{:2d} {:2d} {:2d} {:2d} {:16.12f}\n'\
                                    .format(i+1, k+1, l+1, j+1, \
                                    v1234.real/2))


def twrite_u_matrix_cmplx(v2e, fname, shreshold=1.e-10):
    '''write u-matrix in text form of i, k, l, j, v2e_(ij)(kl).
    '''
    with open(fname, 'w') as f:
        for i, v1 in enumerate(v2e):
            for j, v12 in enumerate(v1):
                for k, v123 in enumerate(v12):
                    if i == k:
                        continue
                    for l, v1234 in enumerate(v123):
                        if j == l:
                            continue
                        if numpy.abs(v1234) > shreshold:
                            # order of cd_i cd_k c_l c_j
                            # absorbing 1/2 factor
                            f.write(('{:2d} {:2d} {:2d} {:2d}'+
                                    ' {:16.12f} {:16.12f}\n')\
                                    .format(i+1, k+1, l+1, j+1, \
                                    v1234.real/2, v1234.imag/2))


def get_u_sab2cshupdn(imp=1):
    '''get unitary transformation from symmetry-adapted basis
    to complex spherical harmonics basis with
    faster spin index ordered as up and down.
    '''
    # default basis to symmetry adapted basis transformation
    with h5py.File('GPARAM.h5', 'r') as f:
        u_db2sab = f['/IMPURITY_{}/DB_TO_SAB'.format(imp)][()].T

    u_sf2of = get_utrans_orbfast_supdn_to_spin_fast_supdn( \
            u_db2sab.shape[0]).T
    u_sf2sab = u_sf2of.dot(u_db2sab)
    u_sab2sf = u_sf2sab.T.conj()
    return u_sab2sf


def h5gen_embedh_spin_updn(imp=1, lwrt_csh=False, lv2e=False):
    # embedding hamiltonian parameters
    # upper case path indicates fortran convention.
    with h5py.File('EMBED_HAMIL_{}.h5'.format(imp), 'r') as f:
        daalpha = f['/D'][()].T
        lambdac = f['/LAMBDA'][()].T
        h1e     = f['/H1E'][()].T
        if lv2e:
            v2e = f['/V2E'][()].T
        else:
            v2e = None

    u_sab2sf = get_u_sab2cshupdn(imp=imp)
    # basis transformation
    daalpha_ = u_sab2sf.T.dot(daalpha).dot(u_sab2sf.conj())
    lambdac_ = u_sab2sf.T.conj().dot(lambdac).dot(u_sab2sf)
    h1e_ = u_sab2sf.T.conj().dot(h1e).dot(u_sab2sf)

    if v2e is not None:
        v2e_ = unitary_transform_coulomb_matrix(v2e, u_sab2sf)
        # double check
        # v2e_ = numpy.einsum('ijkl,pi,jq,rk,ls->pqrs', v2e, \
        #         u_sab2sf.T.conj(), u_sab2sf, \
        #         u_sab2sf.T.conj(), u_sab2sf)
    else:
        v2e_= None

    if lwrt_csh:
        shutil.copy('EMBED_HAMIL_{}.h5'.format(imp), \
                'EMBED_HAMIL_CSH_{}.h5'.format(imp))
        with h5py.File('EMBED_HAMIL_CSH_{}.h5'.format(imp), 'a') as f:
            f['/D'][()] = daalpha_.T
            f['/LAMBDA'][()] = lambdac_.T
            f['/H1E'][()] = h1e_.T
            if v2e_ is not None:
                f['/V2E'][()] = v2e_.T
    return h1e_, lambdac_, daalpha_, v2e_


def h5wrt_rembed_hamil(h1e, lambdac, daalpha, v2e, imp=1):
    '''write the real version of the embedding Hamiltonian in hdf5 format.
    '''
    h1e, lambdac, daalpha, v2e = [get_real_with_chk(a) for a in [h1e, \
            lambdac, daalpha, v2e]]
    fdst = "EMBED_HAMIL_{}r.h5".format(imp)
    if os.path.isfile(fdst):
        with h5py.File(fdst, "a") as f:
            f['/D'][()] = daalpha.T
            f['/LAMBDA'][()] = lambdac.T
            f['/H1E'][()] = h1e.T
    else:
        shutil.copy('EMBED_HAMIL_{}.h5'.format(imp), fdst)
        with h5py.File(fdst, "a") as f:
            del f['/D'], f['/LAMBDA'], f['/H1E'], f['/V2E']
            f['/D'] = daalpha.T
            f['/LAMBDA'] = lambdac.T
            f['/H1E'] = h1e.T
            f['/V2E'] = v2e.T


def get_whole_h1e(h1e, lambdac, daalpha):
    '''get the total one-body part of the embedding hamiltonian
    given the components of impurity, bath and hybridization.
    '''
    norb = h1e.shape[0]
    v1e = block_diag(h1e, -lambdac)
    v1e[:norb, norb:] = daalpha.T
    v1e[norb:, :norb] = daalpha.conj()
    return v1e


def wrt_text_rembed(h1e, lambdac, daalpha, v2e, imp=1, shreshold=1.e-10):
    '''write the real embedding hamiltonian parameters in text format.
    '''
    v1e = get_whole_h1e(h1e, lambdac, daalpha)

    # taking care of f_b fd_a = fd_a f_b + \delta_a,b.
    # v1e += numpy.eye(v1e.shape[0])*lambdac.trace()/lambdac.shape[0]
    # but here the convention is for a typical molecule.
    v1e, v2e = [get_real_with_chk(a) for a in [v1e, v2e]]

    # write one-body part
    twrite_2d_array_real(v1e, 'H1E_{}.INP'.format(imp), shreshold=shreshold)

    # write two-body part
    if v2e is not None:
        twrite_u_matrix_real(v2e, 'V2E_{}.INP'.format(imp),
                shreshold=shreshold)


def wrt_text_cembed(h1e, lambdac, daalpha, v2e, imp=1):
    '''write the complex embedding hamiltonian parameters in text format.
    '''
    v1e = get_whole_h1e(h1e, lambdac, daalpha)

    # taking care of f_b fd_a = fd_a f_b + \delta_a,b.
    # v1e += numpy.eye(v1e.shape[0])*lambdac.trace()/lambdac.shape[0]
    # but here the convention is for a typical molecule.

    # write one-body part
    twrite_2d_array_cmplx(v1e, 'H1E_{}.INP'.format(imp), shreshold=1.e-10)

    # write two-body part
    if v2e is not None:
        twrite_u_matrix_cmplx(v2e, 'V2E_{}.INP'.format(imp))


def get_dm_cshupdn(imp=1):
    '''get the density matrix in complex spherical harmonics basis
    with faster spin index ordered as up, down.
    '''
    with h5py.File('EMBED_HAMIL_RES_{}.h5'.format(imp), 'r') as f:
        dm = f['/DM'][()].T
    u_sab2sf = get_u_sab2cshupdn(imp=imp)
    u2_sab2sf = block_diag(u_sab2sf, u_sab2sf)
    dm = u2_sab2sf.T.dot(dm).dot(u2_sab2sf.conj())
    if numpy.max(numpy.abs(dm.imag)) < 1.e-10:
        dm = dm.real
    numpy.savetxt('dm.dat', dm)


def get_real_with_chk(a, threshold=1.e-10):
    '''return real version of a with a warning if the imaginary part
    is larger than the threshold.
    '''
    if a is None:
        return a
    max_imag = numpy.max(numpy.abs(a.imag))
    if max_imag > 1.e-10:
        print(' maximal imaginary part = {}'.format(max_imag))
    return a.real


def get_res_idmrg(imp=1):
    '''get density matrix and energy from idmrg calculation.
    '''
    # read in density matrix (real)
    fname = "GDMRG_{}.OUT".format(imp)
    dm = numpy.loadtxt(fname)
    if dm.shape[1] == dm.shape[0]*2:
        dm = dm.view(complex)

    # get transformation matrix to symmetry adapted basis
    line = open(fname, "r").readline()
    e_mol = float(line.split()[1])

    return dm, e_mol


def get_res_rspci_mott_onfly(imp=1):
    '''get density matrix and energy from rspci_mott_onfly calculation.
    '''
    with h5py.File("EMBED_HAMIL_RES_{}r.h5".format(imp), "r") as f:
        dm = f["/DM"][()].T
        e_mol = f["/emol"][0]
    return dm, e_mol


def h5wrt_dm_sab_cmplx(dm, e_mol, imp=1):
    '''transform the density matrix (real) into symmetry adapted basis
    and write in EMBEB_HAMIL_{IMP}.h5 file in complex format.
    '''
    # get the comp_sph_harm to sab transformation
    u_sab2sf = get_u_sab2cshupdn(imp=imp)

    # include the transformation for the bath
    u2_sab2sf = block_diag(u_sab2sf, u_sab2sf)

    # unitary transformation
    dm = u2_sab2sf.conj().dot(dm).dot(u2_sab2sf.T)
    with h5py.File("EMBED_HAMIL_RES_{}.h5".format(imp), 'w') as f:
        f["/DM"] = dm.T
        f["/emol"] = [e_mol]
        f["/dimv"] = [0]


def h5wrt_dm_sab_rc(dm, e_mol, imp=1):
    '''transform the density matrix (real) into symmetry adapted basis
    and write in EMBEB_HAMIL_{IMP}.h5 file in complex format.
    '''
    # cast to complex
    dm = numpy.array(dm, dtype=complex)
    h5wrt_dm_sab_cmplx(dm, e_mol, imp=imp)


def chk_consistent_dm(imp=1):
    '''consistent check for dm with unitary transformations.
    '''
    dm = numpy.loadtxt('dm.dat')
    with h5py.File('EMBED_HAMIL_CSH_{}.h5'.format(imp), 'r') as f:
        daalpha = f['/D'][()].T
        lambdac = f['/LAMBDA'][()].T
        h1e     = f['/H1E'][()].T
    v1e = get_whole_h1e(h1e, lambdac, daalpha)
    res1 = numpy.sum(dm*v1e)

    with h5py.File('EMBED_HAMIL_RES_{}.h5'.format(imp), 'r') as f:
        dm = f['/DM'][()].T
    with h5py.File('EMBED_HAMIL_{}.h5'.format(imp), 'r') as f:
        daalpha = f['/D'][()].T
        lambdac = f['/LAMBDA'][()].T
        h1e     = f['/H1E'][()].T
    v1e = get_whole_h1e(h1e, lambdac, daalpha)
    res2 = numpy.sum(dm*v1e)
    print(res1, res2)


def embed_hamil_to_ctext(imp=1):
    with h5py.File("EMBED_HAMIL_{}.h5".format(imp), 'r') as f:
        daalpha = f['/D'][()].T
        lambdac = f['/LAMBDA'][()].T
        h1e     = f['/H1E'][()].T
        v2e     = f['/V2E'][()].T
    wrt_text_cembed(h1e, lambdac, daalpha, v2e, imp=imp)



if __name__ == '__main__':
    _ = h5gen_embedh_spin_updn(lwrt_csh=True, lv2e=True)
