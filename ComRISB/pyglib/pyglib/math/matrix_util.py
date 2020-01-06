"""
Some commonly-used functions of matrix.
"""
import numpy as np
from numpy import linalg as LA
import scipy.sparse as sps
from scipy.sparse import csr_matrix,csc_matrix
from itertools import product

def get_random_hermitian_matrix(n):
    '''
    Generate random hermitian matrix(n, n)
    '''
    x = np.random.rand(n, n) + 1.j * np.random.rand(n, n)
    return x + np.conj(x.T)


def get_random_unitary_matrix(n):
    '''
    Generate random unitary matrix.
    '''
    x = get_random_hermitian_matrix(n)
    w, v = np.linalg.eigh(x)
    return v


def get_matrices_bk_list(matrices, rtol=1.e-12):
    '''
    Get the block dimensions of a group of matrices.
    '''
    l = matrices[0].shape[0]
    l_list = []
    lbase = 0
    for i in range(l):
        is_bk = True
        if i < l - 1:
            for matrix in matrices:
                if np.amax(np.abs(matrix[lbase : i + 1, i + 1 :])) > rtol:
                    is_bk = False
                    break
                if np.amax(np.abs(matrix[i + 1 :, lbase : i + 1])) > rtol:
                    is_bk = False
                    break
        if is_bk:
            l_list.append(i + 1 - lbase)
            lbase = i + 1
    return l_list


def trans_orbital_fast_to_spin_fast(mat, lback=False):
    '''
    Convert a matrix from orbital-fast basis to spin-fast basis or backward
    if lback = True.
    '''
    dim = len(mat)
    U = np.zeros((dim, dim), dtype=int)  # <orbital-fast | spin fast>
    i = 0
    for j in range(0, dim, 2):
        U[i][j] = 1
        i += 1
    for j in range(1, dim, 2):
        U[i][j] = 1
        i += 1
    if lback:
        return np.dot(U, np.dot(mat, U.T))
    else:
        return np.dot(U.T, np.dot(mat, U))


def Rtrans_spin_block(mat):
    '''
    Convert a matrix from spin-down-spin-up block to spin-up-spin-down block
    or backward.
    '''
    dim = len(mat)
    U = np.zeros((dim, dim), dtype=int)
    i = 0
    for j in range(dim / 2, dim):
        U[i][j] = 1
        i += 1
    for j in range(0, dim / 2):
        U[i][j] = 1
        i += 1
    return np.dot(mat, U)


def get_random_unitary_matrix_spin_deg(n, spin_fast=True):
    '''
    Generate spin-degenerate random unitary matrix in the spin-faster basis.
    '''
    u = get_random_unitary_matrix(n)
    from scipy.linalg import block_diag
    u2 = block_diag(u, u)
    if spin_fast:
        u2 = trans_orbital_fast_to_spin_fast(u2)
    return u2


def set_eigen_space(vals):
    '''
    Identify the eigen-space using eigen-values vals.
    '''
    sec_id = [0]
    num_vals = len(vals)
    for i in range(1, num_vals):
        if abs((vals[i] - vals[i - 1]) / max(vals[i - 1], 1.e-10)) > 1.e-4:
            sec_id.append(i)
    sec_id.append(num_vals)
    return sec_id


def sym_array(X, R_list):
    '''
    Symmetrize matrix X.
    '''
    Ham = np.zeros_like(X)
    for i, R in enumerate(R_list):
        if isinstance(R, np.ndarray):
            Ham += R.dot(X).dot(R.conj().T)
        elif isinstance(R, csr_matrix):
            Ham += R.dot(X).dot(R.todense().conj().T)
        else:
            raise AssertionError("Unsupported data type of {}!".\
                    format(R.__class__))
    return Ham / len(R_list)


def sym_csr_matrix(X, R_list):
    '''
    Symmetrize csr_matrix X.
    '''
    from scipy.sparse import csr_matrix
    X_sym = csr_matrix(X.shape, dtype=X.dtype)
    for i, R in enumerate(R_list):
        X_sym += R * X * R.getH()
    return X_sym / len(R_list)


def shrink_label(labels, sec_id, method="average"):
    '''
    Shrink labels for eigen-space.
    method = 'average' ot 'sum'.
    '''
    if labels is None:
        return None
    s_labels = []
    for i in range(len(sec_id) - 1):
        if sec_id[i + 1] <= sec_id[i]:
            continue
        tmp = np.sum(labels[sec_id[i]:sec_id[i + 1]])
        if method == "average":
            tmp /= sec_id[i + 1] - sec_id[i]
        s_labels.append(tmp)
    return np.array(s_labels)


def Rtrans_U_supdn_to_spin_fast_sdnup(U):
    '''
    U = U * untary_{orbital fast (up, down)} {spin fast (down, up)}
    '''
    dim = len(U)
    V = np.zeros((dim, dim), dtype=int)
    i = 0
    for j in range(1, dim, 2):
        V[i][j] = 1
        i += 1
    for j in range(0, dim, 2):
        V[i][j] = 1
        i += 1
    return np.dot(U, V)


def get_utrans_orbfast_supdn_to_spin_fast_supdn(n):
    '''
    get a unitary transformation from orbital fast spin up-down
    basis to spin fast spin up-down basis.
    '''
    v = np.zeros((n, n), dtype=int)
    i = 0

    # spin up block
    for j in range(0, n, 2):
        v[i][j] = 1
        i += 1

    # spin down block
    for j in range(1, n, 2):
        v[i][j] = 1
        i += 1
    return v


def trans_JJ_to_CH_sup_sdn(L):
    '''
    trafoso provides transformation matrices from
    |L,1/2,mL,mS> (L=0,1,2,3, mS=-1/2,1/2) basis to
    basis |J,L,S,mJ>, J=L-1/2, L+1/2
    H. Watanabe 'Operator Methods in Ligand Field Theory'
    Prentice Hall, 1966, Table 1.8-1.
    ordering because of the convention used in WIEN is:
                       mS=1/2        mS=-1/2
                     -L .......L  -L ...... L     (2*(2L+1) columns)
            -(L-1/2)
               .
    J=L-1/2    .
               .
             (L-1/2)
             -L-1/2
               .
    J=L+1/2    .
               .
              L+1/2
    '''
    cf = np.zeros((2 * (2 * L + 1), 2 * (2 * L + 1)))
    if L == 0:
        cf[0, 1] = 1.0
        cf[1, 0] = 1.0
    else:
        k1 = -1
        for ms in range(-1, 2, 2):
            ams = -ms / 2.
            for ml in range(-L, L + 1):
                k1 = k1 + 1
                k2 = -1
                for mj in range(-2 * L + 1, 2 * L, 2):  # L-1/2 states
                    amj = mj / 2.
                    k2 = k2 + 1
                    d = amj - ml - ams
                    if abs(d) < 0.0001:
                        if ms == 1:
                            cf[k2, k1] = - \
                                np.sqrt((L + 0.5 + amj) / (2 * L + 1))
                        else:
                            cf[k2, k1] = np.sqrt((L + 0.5 - amj) / (2 * L + 1))
                for mj in range(-2 * L - 1, 2 * L + 2, 2):  # L+1/2 states
                    amj = mj / 2.
                    k2 = k2 + 1
                    d = amj - ml - ams
                    if abs(d) < 0.0001:
                        if ms == 1:
                            cf[k2, k1] = np.sqrt((L + 0.5 - amj) / (2 * L + 1))
                        else:
                            cf[k2, k1] = np.sqrt((L + 0.5 + amj) / (2 * L + 1))
    return cf


def matrix_basis_expand(mat, basis):
    '''
    Expand mat in basis.
    '''
    res = np.zeros_like(mat)
    for b in basis:
        res += np.vdot(b, mat) * b
    return res


def get_loewner_matrix_ev(e_vals, f, fp):
    '''
    Loewner matrix, ref: positive definite matrices, Bhatia, p.154.
    f: function; fp: first derivative.
    '''
    loewner = []
    for i, ei in enumerate(e_vals):
        _loewner = []
        for j, ej in enumerate(e_vals):
            if np.abs(ei - ej) < 1.e-10:
                res = fp(ei)
            else:
                res = (f(ei) - f(ej))/(ei - ej)
            _loewner.append(res)
        loewner.append(_loewner)
    return np.array(loewner)


def yield_derivative_f_matrix(a, h_list, f, fp):
    '''
    a = \sum_{c} {c*h_list}. fp(x) = \par f(x) / \par x
    Yields \par f(a) / \par c.
    '''
    w, v = LA.eigh(a)
    vherm = np.conj(v.T)
    loewner = get_loewner_matrix_ev(w, f, fp)
    for h in h_list:
        _h = np.dot(vherm, np.dot(h, v))
        yield np.dot(v, np.dot(loewner*_h, vherm))


def get_func_matrix_eigh(A, func, alpha=1.0):
    '''
    Compute matrix function of alpha*A.

    Parameters
    ----------
    A : (N, N) array_like or sparse matrix
        Matrix to be exponentiated.
    alpha : float or complex, optional
        Scale factor. Default is 1.0
    func : string.
        Type of function. Should be one of

            - 'exp' - exponential function.
            - 'log' - natural logarithm.
            - 'sqrt' - square root.
            - 'abs' - absolute value function.

        Otherwise UndefinedFunction exception will be raised.

    Returns
    -------
    res : (N, N) ndarray
        Matrix function of `alpha*A`.
    '''
    if A.shape == (1, 1):
        vals, vecs = np.array([A[0,0]]), np.array([[1.0]], dtype=np.complex)
    else:
        if isinstance(A, csr_matrix) or isinstance(A, csc_matrix):
            A = A.todense()
        vals, vecs = np.linalg.eigh(A)
    if func == 'exp':
        res = np.dot(np.dot(vecs,np.diag(np.exp(alpha * vals))),
                np.conj(vecs).T)
    elif func == 'log':
        # regularize if needed.
        _min = np.min(vals)
        if _min < 1.e-20:
            vals += 1.e-20-_min
        res = np.dot(np.dot(vecs,np.diag(np.log(alpha * vals))),
                np.conj(vecs).T)
    elif func == 'sqrt':
        # regularize if needed.
        _min = np.min(vals)
        if _min < 1.e-20:
            vals += 1.e-20-_min
        res = np.dot(np.dot(vecs,np.diag(np.sqrt(alpha * vals))),
                np.conj(vecs).T)
    elif func == 'abs':
        res = np.dot(np.dot(vecs,np.diag(np.abs(alpha * vals))),
                np.conj(vecs).T)
    else:
        raise Exception('UndefinedFunction')
    return res


def get_matrix_trace(A):
    if A.shape == (1, 1):
        tr = A[0,0]
    elif A.shape == (1,):
        tr = A[0]
    else:
        if sps.issparse(A):
            tr = np.sum(A.diagonal())
        else:
            tr = np.trace(A)
    return tr


def get_loewdin_orthnorm(a):
    '''Get the Loewdin symmetry-adapted orthonormalization via
    singular-value decomposition.
    '''
    U, s, V = np.linalg.svd(a)
    return U.dot(V)


def unitary_transform_coulomb_matrix(a, u):
    '''Perform a unitary transformation (u) on the Coulomb matrix (a).
    '''
    a_ = np.asarray(a).copy()
    m_range = range(a.shape[0])
    for i,j in product(m_range, m_range):
        a_[i,j,:,:] = u.T.conj().dot(a_[i,j,:,:].dot(u))
    a_ = a_.swapaxes(0,2).swapaxes(1,3)
    for i,j in product(m_range, m_range):
        a_[i,j,:,:] = u.T.conj().dot(a_[i,j,:,:].dot(u))
    return a_


def get_hamilt_matrix_from_ev(w, v):
    '''given v_{alpha, n} and w_{n}, get back the Hamiltonian matrix.
    Here alpha is the basis orbital index and n the band index.
    '''
    vh = v.T.conj()
    for i in range(v.shape[0]):
        v[i] *= w
    return v.dot(vh)
