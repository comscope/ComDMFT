from __future__ import print_function
# modified based on U_matrix.py from TRIQS.

from itertools import product
import numpy as np
from scipy.misc import factorial as fact
from pyglib.math.matrix_util import unitary_transform_coulomb_matrix


def U_matrix(mode, l, radial_integrals=None, U_int=None, J_hund=None, T=None):
    if 'slater' in mode:
        U_matrix = U_matrix_slater(l, radial_integrals, U_int, J_hund)
    elif 'kanamori' in mode:
        U_matrix = U_matrix_kanamori(2*l+1, U_int, J_hund)
    else:
        raise NameError(" unsupported mode!")
    u_avg, j_avg = get_average_uj(U_matrix)

    # add spin-components
    norb = U_matrix.shape[0]
    norb2 = norb*2
    Ufull_matrix = np.zeros((norb2,norb2,norb2,norb2), dtype=np.complex)
    if T is not None:
        # spin block
        Ufull_matrix[:norb,:norb,:norb,:norb] = U_matrix
        Ufull_matrix[norb:,norb:,norb:,norb:] = U_matrix
        Ufull_matrix[:norb,:norb,norb:,norb:] = U_matrix
        Ufull_matrix[norb:,norb:,:norb,:norb] = U_matrix
        print(" u-matrix: nnz in compl_sph_harm = {}".format(
                np.count_nonzero(np.abs(Ufull_matrix)>1.e-10)))
        Ufull_matrix = unitary_transform_coulomb_matrix(Ufull_matrix, T)
    else: # spin index fast
        Ufull_matrix[::2,::2,::2,::2] = U_matrix  # up, up
        Ufull_matrix[1::2,1::2,1::2,1::2] = U_matrix # dn, dn
        Ufull_matrix[::2,::2,1::2,1::2] = U_matrix # up, dn
        Ufull_matrix[1::2,1::2,::2,::2] = U_matrix # dn, up

    print(" u-matrix: nnz in final basis = {}".format(
            np.count_nonzero(np.abs(Ufull_matrix)>1.e-10)))
    return Ufull_matrix, u_avg, j_avg


# The interaction matrix in desired basis.
# U^{spherical}_{m1 m4 m2 m3} =
# \sum_{k=0}^{2l} F_k angular_matrix_element(l, k, m1, m2, m3, m4)
# H = \frac{1}{2} \sum_{ijkl,\sigma \sigma'} U_{ikjl}
# a_{i \sigma}^\dagger a_{j \sigma'}^\dagger a_{l \sigma'} a_{k \sigma}.
def U_matrix_slater(l, radial_integrals=None, U_int=None, J_hund=None):
    r"""
    Calculate the full four-index U matrix being given either
    radial_integrals or U_int and J_hund.
    The convetion for the U matrix is that used to construct
    the Hamiltonians, namely:

    .. math:: H = \frac{1}{2} \sum_{ijkl,\sigma \sigma'} U_{ikjl}
            a_{i \sigma}^\dagger a_{j \sigma'}^\dagger
            a_{l \sigma'} a_{k \sigma}.

    Parameters
    ----------
    l : integer
        Angular momentum of shell being treated
        (l=2 for d shell, l=3 for f shell).
    radial_integrals : list, optional
                       Slater integrals [F0,F2,F4,..].
                       Must be provided if U_int and J_hund are not given.
                       Preferentially used to compute the U_matrix
                       if provided alongside U_int and J_hund.
    U_int : scalar, optional
            Value of the screened Hubbard interaction.
            Must be provided if radial_integrals are not given.
    J_hund : scalar, optional
             Value of the Hund's coupling.
             Must be provided if radial_integrals are not given.

    Returns
    -------
    U_matrix : float numpy array
               The four-index interaction matrix in the chosen basis.
    """

    # Check all necessary information is present and consistent
    if radial_integrals is None and (U_int is None and J_hund is None):
        raise ValueError("U_matrix: provide either the radial_integrals" +
                " or U_int and J_hund.")
    if radial_integrals is None and (U_int is not None and J_hund is not None):
        radial_integrals = U_J_to_radial_integrals(l, U_int, J_hund)
    if radial_integrals is not None and \
            (U_int is not None and J_hund is not None):
        if len(radial_integrals)-1 != l:
            raise ValueError("U_matrix: inconsistency in l" +
                    " and number of radial_integrals provided.")
        if not np.allclose(radial_integrals,
                U_J_to_radial_integrals(l, U_int, J_hund)):
            print(" Warning: U_matrix: radial_integrals provided\n"+
            " do not match U_int and J_hund.\n"+
            " Using radial_integrals to calculate U_matrix.")

    # Full interaction matrix
    # Basis of spherical harmonics Y_{-2}, Y_{-1}, Y_{0}, Y_{1}, Y_{2}
    # U^{spherical}_{m1 m4 m2 m3} = \sum_{k=0}^{2l} F_k
    # angular_matrix_element(l, k, m1, m2, m3, m4)
    U_matrix = np.zeros((2*l+1,2*l+1,2*l+1,2*l+1),dtype=np.float)

    m_range = range(-l,l+1)
    for n, F in enumerate(radial_integrals):
        k = 2*n
        for m1, m2, m3, m4 in product(m_range,m_range,m_range,m_range):
            U_matrix[m1+l,m3+l,m2+l,m4+l] += \
                    F * angular_matrix_element(l,k,m1,m2,m3,m4)
    return U_matrix


def U_matrix_kanamori(n_orb, U_int, J_hund):
    r"""
    Calculate the Kanamori U and Uprime matrices.

    Parameters
    ----------
    n_orb : integer
            Number of orbitals in basis.
    U_int : scalar
            Value of the screened Hubbard interaction.
    J_hund : scalar
             Value of the Hund's coupling.

    Returns
    -------
    U_matrix : float numpy array
               The four-index interaction matrix in the chosen basis.
    """

    U_matrix  = np.zeros((n_orb,n_orb,n_orb,n_orb),dtype=np.float)
    m_range = range(n_orb)
    for m,mp in product(m_range,m_range):
        if m == mp:
            U_matrix[m,m,mp,mp] = U_int
        else:
            U_matrix[m,m,mp,mp] = U_int - 2.0*J_hund
            U_matrix[m,mp,mp,m] = J_hund
            U_matrix[m,mp,m,mp] = J_hund
    return U_matrix


def U_J_to_radial_integrals4(l,U,J):
    label_to_l = {'s':0, 'p':1, 'd':2, 'f':3}
    if isinstance(l,basestring):
        l = label_to_l[l_label]
    f_list = U_J_to_radial_integrals(l, U, J)
    for i in range(len(f_list),4):
        f_list.append(0.)
    return f_list


# Convert U,J -> radial integrals F_k
def U_J_to_radial_integrals(l, U_int, J_hund):
    r"""
    Determine the radial integrals F_k from U_int and J_hund.

    Parameters
    ----------
    l : integer
        Angular momentum of shell being treated
        (l=2 for d shell, l=3 for f shell).
    U_int : scalar
            Value of the screened Hubbard interaction.
    J_hund : scalar
             Value of the Hund's coupling.

    Returns
    -------
    radial_integrals : list
                       Slater integrals [F0,F2,F4,..].
    """

    F = np.zeros((l+1),dtype=np.float)
    F[0] = U_int
    if l == 0:
        pass
    elif l == 1:
        F[1] = J_hund * 5.0
    elif l == 2:
        F[1] = J_hund * 14.0 / (1.0 + 0.625)
        F[2] = 0.625 * F[1]
    elif l == 3:
        F[1] = 6435.0 * J_hund / (286.0 + 195.0 * 0.668 + 250.0 * 0.494)
        F[2] = 0.668 * F[1]
        F[3] = 0.494 * F[1]
    else:
        raise ValueError(
            " U_J_to_radial_integrals: implemented only for l=0,1,2,3")
    return F


# Convert radial integrals F_k -> U,J
def radial_integrals_to_U_J(l, F):
    r"""
    Determine U_int and J_hund from the radial integrals.

    Parameters
    ----------
    l : integer
        Angular momentum of shell being treated
        (l=2 for d shell, l=3 for f shell).
    F : list
        Slater integrals [F0,F2,F4,..].

    Returns
    -------
    U_int : scalar
            Value of the screened Hubbard interaction.
    J_hund : scalar
             Value of the Hund's coupling.

    """

    U_int = F[0]
    if l == 0:
        J_Hund = 0.
    elif l == 1:
        J_hund = F[1]/5.0
    elif l == 2:
        J_hund = F[1]*(1.0 + 0.625)/14.0
    elif l == 3:
        J_hund = F[1] * (286.0 + 195.0*0.668 + 250.0*0.494)/6435.0
    else: raise ValueError("radial_integrals_to_U_J:"+
            " implemented only for l=2,3")

    return U_int,J_hund


# Angular matrix elements of particle-particle interaction
# (2l+1)^2 ((l 0) (k 0) (l 0))^2 \sum_{q=-k}^{k} (-1)^{m1+m2+q}
# ((l -m1) (k q) (l m3)) ((l -m2) (k -q) (l m4))
def angular_matrix_element(l, k, m1, m2, m3, m4):
    r"""
    Calculate the angular matrix element

    .. math::
       (2l+1)^2
       \begin{pmatrix}
            l & k & l \\
            0 & 0 & 0
       \end{pmatrix}^2
       \sum_{q=-k}^k (-1)^{m_1+m_2+q}
       \begin{pmatrix}
            l & k & l \\
         -m_1 & q & m_3
       \end{pmatrix}
       \begin{pmatrix}
            l & k  & l \\
         -m_2 & -q & m_4
       \end{pmatrix}.

    Parameters
    ----------
    l : integer
    k : integer
    m1 : integer
    m2 : integer
    m3 : integer
    m4 : integer

    Returns
    -------
    ang_mat_ele : scalar
                  Angular matrix element.

    """
    ang_mat_ele = 0
    for q in range(-k,k+1):
        ang_mat_ele += three_j_symbol((l,-m1),(k,q),(l,m3))* \
                three_j_symbol((l,-m2),(k,-q),(l,m4))* \
                (-1.0 if (m1+q+m2) % 2 else 1.0)
    ang_mat_ele *= (2*l+1)**2 * (three_j_symbol((l,0),(k,0),(l,0))**2)
    return ang_mat_ele


# Wigner 3-j symbols
# ((j1 m1) (j2 m2) (j3 m3))
def three_j_symbol(jm1, jm2, jm3):
    r"""
    Calculate the three-j symbol
    .. math::
       \begin{pmatrix}
        l_1 & l_2 & l_3\\
        m_1 & m_2 & m_3
       \end{pmatrix}.
    Parameters
    ----------
    jm1 : tuple of integers
          (j_1 m_1)
    jm2 : tuple of integers
          (j_2 m_2)
    jm3 : tuple of integers
          (j_3 m_3)
    Returns
    -------
    three_j_sym : scalar
                  Three-j symbol.
    """
    j1, m1 = jm1
    j2, m2 = jm2
    j3, m3 = jm3

    if (m1+m2+m3 != 0 or
        m1 < -j1 or m1 > j1 or
        m2 < -j2 or m2 > j2 or
        m3 < -j3 or m3 > j3 or
        j3 > j1 + j2 or
        j3 < abs(j1-j2)):
        return .0

    three_j_sym = -1.0 if (j1-j2-m3) % 2 else 1.0
    three_j_sym *= np.sqrt(fact(j1+j2-j3)*fact(j1-j2+j3)* \
            fact(-j1+j2+j3)/fact(j1+j2+j3+1))
    three_j_sym *= np.sqrt(fact(j1-m1)*fact(j1+m1)*fact(j2-m2)* \
            fact(j2+m2)*fact(j3-m3)*fact(j3+m3))

    t_min = max(j2-j3-m1,j1-j3+m2,0)
    t_max = min(j1-m1,j2+m2,j1+j2-j3)

    t_sum = 0
    for t in range(t_min,t_max+1):
        t_sum += (-1.0 if t % 2 else 1.0)/(fact(t)*fact(j3-j2+m1+t)* \
                fact(j3-j1-m2+t)*fact(j1+j2-j3-t)*fact(j1-m1-t)*fact(j2+m2-t))

    three_j_sym *= t_sum
    return three_j_sym


def get_average_uj(v2e):
    m_range = range(v2e.shape[0])
    u_avg = 0; j_avg = 0
    isum_u = 0; isum_j = 0
    for i, j in product(m_range, m_range):
        u_avg += v2e[i, i, j, j]
        isum_u += 1
        if i != j:
            j_avg += v2e[i, i, j, j] - v2e[i, j, j, i]
            isum_j += 1
    u_avg /= isum_u
    if isum_j > 0:
        j_avg = u_avg - j_avg/isum_j
    return u_avg, j_avg


def get_v2e_list(lhub, l_list, imap_list, utrans_list, u_list=None,
        j_list=None, f_list=None):
    mode_list = ['manual', 'slater-condon', 'kanamori', 'slater-condon']
    if lhub > 0:
        v2e_list = []
        u_avg_list = []
        j_avg_list = []

        for i, imap in enumerate(imap_list):
            if i > imap:
                v2e_list.append(v2e_list[imap])
                u_avg_list.append(u_avg_list[imap])
                j_avg_list.append(j_avg_list[imap])
                continue
            assert len(l_list[i]) == 1, " more than one l with lhub>0!"
            l_imp = l_list[i][0]
            utrans = utrans_list[i]

            if lhub in [1, 2]:
                v2e, u_avg, j_avg = U_matrix(mode_list[lhub], l_imp,
                        U_int=u_list[i],
                        J_hund=j_list[i], T=utrans)
            else:
                v2e, u_avg, j_avg = U_matrix(mode_list[lhub], l_imp,
                        radial_integrals=f_list[i][:l_imp+1], T=utrans)

            v2e_list.append(v2e)
            u_avg_list.append(u_avg)
            j_avg_list.append(j_avg)
    else:
        v2e_list = None
        u_avg_list = None
        j_avg_list = None

    return v2e_list, u_avg_list, j_avg_list


if __name__ == '__main__':
    l = 2
    coulomb_matrix = U_matrix('slater', l, U_int=6.7, J_hund=0.67)
    m_range = range(2*l+1)
    with open('coulomb_matrixp.txt','w') as f:
        for m1,m2,m3,m4 in product(m_range,m_range,m_range,m_range):
            if np.abs(coulomb_matrix[m1,m2,m3,m4]) < 1.e-8:
                continue
            f.write('{} {} {} {}  {}\n'.format(m1+1,m2+1,m3+1,m4+1, \
                    coulomb_matrix[m1,m2,m3,m4]))
