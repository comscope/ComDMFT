# -*- coding: utf-8 -*-
from __future__ import print_function
from pyglib.symm.unitary import comp_sph_harm_to_relativistic_harm

'''
Get the matrix representation of angular momentum operator
in complex spherical harmonics (CSH) basis.
The orbital index is faster and spin goes up first where applicable.
'''

import numpy as np
from scipy.linalg import block_diag
from pyglib.math.matrix_util import trans_orbital_fast_to_spin_fast


def get_J_generator(l_list, iso):
    '''
    Get rotation generator J vector from a l_list with CSH basis.
    For iso=2, the fast index of the basis is orbital,
    and the spin goes up first then down.
    '''
    if iso == 1:
        # need L-vector for one spin-block only
        # because of no spin-orbit interaction.
        return get_L_vector(l_list)
    else:
        return get_J_vector(l_list)


def get_J_vector(l_list):
    '''get the total angular momentum vector in CSH basis,
    with orbital index faster and spin running through [up, down].
    '''
    Lx, Ly, Lz = get_L_vector(l_list, iso=2)
    Sx, Sy, Sz = get_S_vector(l_list)
    return [Lx + Sx, Ly + Sy, Lz + Sz]


def get_JU_relat_sph_harm_cg(l_list):
    '''get J vector in relativisitic spherical harmonics (RSH) basis,
    and the unitary transformation U from complex spherical harmonics
    (CSH) to RSH using clebsch-gordan coefficients.
    '''
    for i, _l in enumerate(l_list):
        _Jx, _Jy, _Jz =  get_J_vector(_l)
        v = comp_sph_harm_to_relativistic_harm((2*_l+1)*2)
        _Jx, _Jy, _Jz = v.T.conj().dot(_Jx).dot(v), \
                v.T.conj().dot(_Jy).dot(v), v.T.conj().dot(_Jz).dot(v)
        # check Jz
        for j in range(2*_l):
            jz = j-(_l-0.5)
            if abs(jz-_Jz[j,j]) > 1.e-6:
                raise ValueError("jz = {} vs expected {}!".format(\
                        _Jz[j], jz))
        for j in range(2*_l+2):
            jz = j-(_l+0.5)
            if abs(jz-_Jz[j+2*_l,j+2*_l]) > 1.e-6:
                raise ValueError("jz = {} vs expected {}!".format(\
                        _Jz[j+2*_l,j+2*_l], jz))
        if i == 0:
            Jx, Jy, Jz = _Jx, _Jy, _Jz
            u_trans = v
        else:
            Jx, Jy, Jz, u_trans = block_diag(Jx, _Jx), block_diag(Jy, _Jy), \
                    block_diag(Jz, _Jz), block_diag(u_trans, v)
    return [Jx, Jy, Jz], u_trans


def get_JU_relat_sph_harm_random_phase(l_list):
    '''get J vector in relativisitic spherical harmonics (RSH) basis,
    and the unitary transformation U from complex spherical harmonics
    (CSH) to RSH.
    There is a random phase associated RSH due to diagonalization.
    '''
    for i, _l in enumerate(l_list):
        _Jx, _Jy, _Jz =  get_J_vector(_l)
        Jsq = _Jx.dot(_Jx) + _Jy.dot(_Jy) + _Jz.dot(_Jz)
        w, v = np.linalg.eigh(Jsq)

        # from j(j+1) to j
        w = map(lambda x: np.sqrt(x+0.25)-0.5, w)

        if _l > 0:
            n1 = int(2*(_l - .5) + 1.1)
            if np.max(np.abs(w[:n1] - w[0])) > 1.e-6 or \
                    np.max(np.abs(w[n1:] - w[n1])) > 1.e-6 or \
                    np.abs(w[n1] - w[0] - 1) > 1.e-6:
                raise_j_list_error(w, 'j_list')

            # l-1/2 block
            Jz1 = v[:,:n1].T.conj().dot(_Jz).dot(v[:,:n1])
            w1, v1 = np.linalg.eigh(Jz1)
            if np.max(np.abs(w1[1:] - w1[:-1] - 1.)) > 1.e-6:
                raise_j_list_error(w1, 'jz1_list')
            v[:, :n1] = v[:, :n1].dot(v1)
        else:
            n1 = 0

        # l+1/2 block
        Jz1 = v[:,n1:].T.conj().dot(_Jz).dot(v[:,n1:])
        w1, v1 = np.linalg.eigh(Jz1)
        if np.max(np.abs(w1[1:] - w1[:-1] - 1.)) > 1.e-6:
            raise_j_list_error(w1, 'jz2_list')
        v[:, n1:] = v[:, n1:].dot(v1)
        _Jx, _Jy, _Jz = v.T.conj().dot(_Jx).dot(v), \
                v.T.conj().dot(_Jy).dot(v), v.T.conj().dot(_Jz).dot(v)
        if i == 0:
            Jx, Jy, Jz = _Jx, _Jy, _Jz
            u_trans = v
        else:
            Jx, Jy, Jz, u_trans = block_diag(Jx, _Jx), block_diag(Jy, _Jy), \
                    block_diag(Jz, _Jz), block_diag(u_trans, v)
    return [Jx, Jy, Jz], u_trans


def raise_j_list_error(j_list, head):
    msg = 'error in {}:\n'.format(head)
    for k, j in enumerate(j_list):
        msg += ' {:6.2f}'.format(j)
        if k % 5 == 4:
            msg += '/n'
    raise ValueError(msg)


def get_L_vector(l, iso=1):
    '''
    Get matrix representation of L-vector.
    '''
    if iso == 1:
        # for spin-up block
        try:
            for i, _l in enumerate(l):
                _Lx, _Ly, _Lz = get_L_vector(_l)
                if i == 0:
                    Lx, Ly, Lz = _Lx.copy(), _Ly.copy(), _Lz.copy()
                else:
                    Lx, Ly, Lz = block_diag(Lx, _Lx), block_diag(Ly, _Ly), \
                        block_diag(Lz, _Lz)
            return [Lx, Ly, Lz]
        # excepct single l
        except TypeError:
            Lz = get_matrix_Lz_CSH(l)
            Lp = get_matrix_Lp_CSH(l)
            _, Lx, Ly = get_other_op(Lz, Lp)
            return [Lx, Ly, Lz]
    else:
        # get spin-up block first
        Lx, Ly, Lz = get_L_vector(l)
        # add spin-down block
        Lx, Ly, Lz = block_diag(Lx, Lx), block_diag(Ly, Ly), \
                block_diag(Lz, Lz)
        return [Lx, Ly, Lz]


def get_Lp_coef(l, m):
    '''
    L^\dagger |l, m> = Lp_coef |l, m+1>.
    '''
    return np.sqrt(l * (l + 1) - m * (m + 1))


def get_diag_Lz_CSH(l):
    '''
    Get diagonal elements of Lz for Complex spherical Harmomnics basis.
    '''
    Lz = []
    for j in range(int(2 * l + 1)):
        Lz.append(-l + j)
    return Lz


def get_matrix_Lz_CSH(l):
    '''
    Get matrix_{z}.
    '''
    return np.diag(get_diag_Lz_CSH(l))


def get_matrix_Lp_CSH(l, m_rising=True):
    '''
    Get matrix L_{\dagger}.
    '''
    lm = int(2 * l + 1.5)
    Lp = np.zeros([lm, lm])
    if m_rising:
        for i in range(int(2 * l + 0.5)):
            Lp[i + 1, i] = get_Lp_coef(l, -l + i)
    else:
        for i in range(int(2 * l + 0.5)):
            Lp[i, i + 1] = get_Lp_coef(l, l - i - 1)
    return Lp


def get_matrix_Sz_CSH_orbital_fast(l_list):
    '''
    Get matrix Sz for Complex spherical Harmomnics basis
    with spin-up + spin-down block.
    '''
    num_lm = np.sum(2*np.asarray(l_list) + 1)

    # spin-up + spin_dn. Wien2k convention.
    Sz = np.diag([ 0.5 for i in range(num_lm)] + \
            [-0.5 for i in range(num_lm)])
    return Sz


def get_matrix_Sp_CSH_spin_fast(l_list):
    '''
    Get matrix Sp for Complex spherical Harmomnics basis with spin-fast
    (up, down) index.
    '''
    num_lm = np.sum(2 * np.asarray(l_list) + 1)
    Sp_sub = get_matrix_Lp_CSH(0.5, m_rising=False)
    Sp_list = [Sp_sub for i in range(num_lm)]
    from scipy.linalg import block_diag
    return block_diag(*Sp_list)


def get_S_vector(l_list):
    '''
    Get S-vector in CSH basis with fast orbital index.
    '''
    Sz = get_matrix_Sz_CSH_orbital_fast(l_list)
    Sp = get_matrix_Sp_CSH_spin_fast(l_list)
    Sp = trans_orbital_fast_to_spin_fast(Sp, lback=True)
    _, Sx, Sy = get_other_op(Sz, Sp)
    return [Sx, Sy, Sz]


def get_other_op(Lz, Lp):
    '''
    Get other related operators given L_z and L^\dagger.
    '''
    Ln = Lp.conj().T
    Lx = (Lp + Ln) / 2.
    Ly = (Lp - Ln) / (2.j)
    return Ln, Lx, Ly


def get_l_list_from_string(code):
    '''
    Gt l list from letters, e.g., "spd" -> [0, 1, 2].
    '''
    l_list = []
    for c in code:
        l = get_l_from_char(c)
        if l >= 0:
            l_list.append(l)
    return l_list


def get_l_from_char(c):
    return {'s': 0, 'p': 1, 'd': 2, 'f': 3}.get(c, -1)


def get_complex_to_real_sph_harm(l):
    '''get the unitary transformation from compex spherical harmonics
    (Condon–Shortley phase convention) to real harmonics
    (https://en.wikipedia.org/wiki/Spherical_harmonics,
    consistent with VASP).
    '''
    mdim = 2*l + 1
    c2r = np.zeros((mdim, mdim), dtype=np.complex)
    for m in range(mdim):
        m_ = m - l
        if m_ > 0:
            c2r[ m_+l, m] = (-1)**m_/np.sqrt(2.)
            c2r[-m_+l, m] = 1./np.sqrt(2.)
        elif m_ == 0:
            c2r[l, l] = 1.
        else:
            c2r[ m_+l, m] = 1.j/np.sqrt(2.)
            c2r[-m_+l, m] = -1.j*(-1)**m_/np.sqrt(2.)
    return c2r


def get_jvec_comp_sph_harm_to_rel_harm(j_csh, l_list):
    '''get the J vector in complex spherical harmonics representation
    to relativistic harmonics representation.
    the transformation matrix is also returned.
    '''
    lmbase = 0
    for il, l in enumerate(l_list):
        lm2 = (l*2+1)*2
        jx = j_csh[0][lmbase:lmbase+lm2, lmbase:lmbase+lm2]
        jy = j_csh[1][lmbase:lmbase+lm2, lmbase:lmbase+lm2]
        jz = j_csh[2][lmbase:lmbase+lm2, lmbase:lmbase+lm2]
        lmbase += lm2
        j2 = jx.dot(jx) + jy.dot(jy) + jz.dot(jz)
        w, v = np.linalg.eigh(j2)
        jx = v.conj().T.dot(jx).dot(v)
        jy = v.conj().T.dot(jy).dot(v)
        jz = v.conj().T.dot(jz).dot(v)
        if il == 0:
            j_rel = [jx, jy, jz]
            u_csh2rel = v
        else:
            j_rel[0] = block_diag(j_rel[0], jx)
            j_rel[1] = block_diag(j_rel[1], jy)
            j_rel[2] = block_diag(j_rel[2], jz)
            u_csh2rel = block_diag(u_csh2rel, v)
    return j_rel, u_csh2rel



if __name__ == "__main__":
    Lx, Ly, Lz = get_L_vector(0)
    print(Lx, Ly, Lz)
