from __future__ import print_function

'''
Tools to analyse the local many-body density matrix (multiplet structure).
'''

try:
    from builtins import range
except:
    pass

import numpy as np
from scipy.linalg import logm
import h5py


def get_rho_histogram(rho, S=None, L=None, J=None, num_ev=0, Rpr_list=None):
    '''
    Get the histogram of the local density matrix with labels.

    This function diagonalizes the reduced local many-body density matrix rho,
    in the order of valence block.
    The label of the resultant eigen-space, such as (averaged) J,
    character chi of the irreducible representation and degeneracy,
    will be computed.
    The final results are ordered accordig to descending order of
    the eigen-values of the eigen-spaces.

    Parameters
    ----------
    rho : 2d array
        reduced local many-body density matrix in sparse matrix format.
    J : csc_matrix
        Total angular momentum :math:`J^{2}` operator in sparse matrix format.
    num_ev : integer
        Number of significant eigen-vectors of rho to be calculated.
    Rpr_list : list.
        List of rotation operations in the local Hilbert space of difference
        valence block.

    Returns
    -------
    vals : array
        Eigen-values of the eigen-spaces of rho.
    j_label : array
        Averaged J values of the eigen-spaces of rho.
    chi_label : array
        Characters of the eigen-spaces of rho as the irreducible
        representation of the rotation group defined by Rpr_list.
    multiplet_degeneracies : array
        Degeneracies of the eigen-spaces of rho.

    '''
    print(" Get rho histogram.")
    if num_ev > rho.shape[0]:
        num_ev = 0
    if num_ev <= 0:
        vals, vecs = np.linalg.eigh(rho)
        vecs = vecs.T
    else:
        from scipy.sparse.linalg import eigsh
        vals, vecs = eigsh(rho, num_ev, which='LM')
        vecs = vecs.T
    # Make sure eigen-values properly ordered.
    idx = vals.argsort()[::-1]
    vals = vals[idx]
    vecs = vecs[idx, :]

    from pyglib.mbody.local_operator_factory import get_label_list
    print(" Get labels.")
    s_label = get_label_list(S, vecs)
    if s_label is not None:
        s_label = np.sqrt(s_label + 0.25) - 0.5

    l_label = get_label_list(L, vecs)
    if l_label is not None:
        l_label = np.sqrt(l_label + 0.25) - 0.5

    j_label = get_label_list(J, vecs)
    if j_label is not None:
        j_label = np.sqrt(j_label + 0.25) - 0.5

    from pyglib.math.matrix_util import set_eigen_space, shrink_label
    idx = set_eigen_space(vals)
    vals = shrink_label(vals, idx, method="sum")
    s_label = shrink_label(s_label, idx, method="average")
    l_label = shrink_label(l_label, idx, method="average")
    j_label = shrink_label(j_label, idx, method="average")
    multiplet_degeneracies = np.array(
        [idx[i + 1] - idx[i] for i in range(len(idx) - 1)])

    if Rpr_list is not None:
        # check commutation
        from pyglib.symm.atom_symm import check_commute_G
        check_commute_G(rho, Rpr_list)

        from pyglib.symm.atom_symm import get_characters_espace, \
                check_sum_chi2_1
        chi_label, _ = get_characters_espace(Rpr_list, vecs)
        check = check_sum_chi2_1(chi_label)
        if check != "OK":
            print(" Warning: chi-sum error!")
    else:
        chi_label = None

    return vals, s_label, l_label, j_label, chi_label, multiplet_degeneracies


def get_ordered_labels(val_list, n_list, s_list, l_list, j_list, \
        chi_list, d_list):
    '''
    Sort according to eigen-values of rho.
    '''

    idx = np.asarray(val_list).argsort()[::-1]
    val_list = np.asarray(val_list)[idx]
    n_list = np.asarray(n_list)[idx]
    if s_list is not None:
        s_list = np.asarray(s_list)[idx]
    if l_list is not None:
        l_list = np.asarray(l_list)[idx]
    if j_list is not None:
        j_list = np.asarray(j_list)[idx]
    if chi_list != []:
        chi_list = np.asarray(chi_list)[idx]
    d_list = np.asarray(d_list)[idx]
    return val_list, n_list, s_list, l_list, j_list, chi_list, d_list


def get_local_histogram(imp, op_list=["S2","L2","J2"], num_ev=0,
        Rpr_list=None, n_bot = None, n_top = None):
    '''
    Get the local histogram with labels.
    '''

    with h5py.File('EMBED_HAMIL_{}.h5'.format(imp), 'r') as f:
        _nbot = f['/nval_bot'][0]
        _ntop = f['/nval_top'][0]
    if n_bot is None:
        n_bot = _nbot
    if n_top is None:
        n_top = _ntop
    n_bot = max(n_bot, _nbot)
    n_top = min(n_top, _ntop)
    print(' n_bot = {}, n_top = {}'.format(n_bot, n_top))

    val_list = []
    s_list = []
    l_list = []
    j_list = []
    n_list = []
    chi_list = []
    d_list = []

    if "S2" in op_list:
        eval_s2 = 0
    else:
        eval_s2 = None
    if "L2" in op_list:
        eval_l2 = 0
    else:
        eval_l2 = None
    if "J2" in op_list:
        eval_j2 = 0
    else:
        eval_j2 = None
    ent_entropy = 0

    from pyglib.mbody.local_operator_factory import get_local_operators
    for ival in range(n_bot, n_top+1):
        print(" Valence = {}".format(ival))
        # Get reduced density matrix.
        with h5py.File('EMBED_HAMIL_ANALYSIS_{}.h5'.format(imp), 'r') as f:
            path = '/valence_block_{}/RHO'.format(ival)
            if path not in f:
                continue
            rho = f[path][...].T

        tr_rho = np.trace(rho)
        print(' trace(rho) = {}'.format(tr_rho))

        # entanglement entropy
        if tr_rho > 1.e-6:
            ent_entropy -= np.trace(rho.dot(logm(rho)))

        print(" Get local operators (J2, ...)")
        op_dict = get_local_operators(imp, ival, op_list)

        if "S2" in op_list:
            S = op_dict["S2"]
            eval_s2 += np.trace(S.dot(rho))
        else:
            S = None

        if "L2" in op_list:
            L = op_dict["L2"]
            eval_l2 += np.trace(L.dot(rho))
        else:
            L = None

        if "J2" in op_list:
            J = op_dict["J2"]
            eval_j2 += np.trace(J.dot(rho))
        else:
            J = None

        if Rpr_list is not None:
            Rpr_list = Rpr_list[n_bot : n_top + 1]

        vals, s_label, l_label, j_label, chi_label, degeneracies = \
                get_rho_histogram(rho, S=S, L=L, J=J, num_ev=num_ev, \
                Rpr_list=Rpr_list)
        val_list.extend(vals)

        d_list.extend(degeneracies)
        if s_label is not None:
            s_list.extend(s_label)
        if l_label is not None:
            l_list.extend(l_label)
        if j_label is not None:
            j_list.extend(j_label)
        if chi_label is not None:
            chi_list.extend(chi_label)
        n_list.extend([ival for i in enumerate(vals)])

    # Convert from <S^2> to <S>
    if eval_s2 is not None:
        eval_s2 = np.sqrt(eval_s2 + 0.25) - 0.5
    if eval_l2 is not None:
        eval_l2 = np.sqrt(eval_l2 + 0.25) - 0.5
    if eval_j2 is not None:
        eval_j2 = np.sqrt(eval_j2 + 0.25) - 0.5

    val_list, n_list, s_list, l_list, j_list, chi_list, d_list = \
            get_ordered_labels(val_list, n_list, s_list, l_list, j_list, \
            chi_list, d_list)

    return val_list, n_list, s_list, l_list, j_list, chi_list, d_list, \
            eval_s2, eval_l2, eval_j2, ent_entropy

