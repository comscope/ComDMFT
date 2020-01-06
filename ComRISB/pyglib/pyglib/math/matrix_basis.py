"""
Given the matrix structure, generate the (Hermitian) matrix basis set.
"""
import numpy as np
from scipy.sparse import csc_matrix


def sigmatomatrixbasis(sigma):
    '''
    Generate Hermitian matrix basis set.
    '''
    matrix_basis = []
    sigma = np.asarray(sigma)
    for element in range(np.max(sigma), 0, -1):
        spots = np.argwhere(sigma == element)
        num_spots = len(spots)
        if num_spots == 0:
            continue
        spots = spots.T
        # Skip if located at lower trigonal block
        if spots[0][0] > spots[1][0]:
            continue
        if spots[0][0] == spots[1][0]:
            value = 1 / np.sqrt(float(num_spots))
            matrix = np.zeros_like(sigma, dtype=np.complex)
            matrix[spots[0], spots[1]] = value
            matrix_basis.append(matrix)
        else:
            # non-zero element
            value = 1 / np.sqrt(float(num_spots * 2))
            matrix = np.zeros_like(sigma, dtype=np.complex)
            matrix[spots[0], spots[1]] = value
            matrix[spots[1], spots[0]] = value
            matrix_basis.append(matrix)
            value = value * 1.j
            matrix = np.zeros_like(sigma, dtype=np.complex)
            matrix[spots[0], spots[1]] = value
            matrix[spots[1], spots[0]] = -value
            matrix_basis.append(matrix)
    return matrix_basis


def listsigmatomatrixbasis(sigma_list):
    matrix_basis_list = []
    for sigma in sigma_list:
        matrix_basis_list.append(SigmaToMatrixBasis(sigma))
    return matrix_basis_list


def matrixstructtobasis(m_struct):
    '''
    Generate general non-Hermitian matrix basis set.
    '''
    matrix_basis = []
    m_struct = np.asarray(m_struct)
    for element in range(np.max(m_struct), 0, -1):
        spots = np.argwhere(m_struct == element)
        num_spots = len(spots)
        if num_spots == 0:
            continue
        spots = spots.T
        value = 1 / np.sqrt(float(num_spots))
        matrix = np.zeros_like(m_struct, dtype=np.complex)
        matrix[spots[0], spots[1]] = value
        matrix_basis.append(matrix)
    return matrix_basis


def listmatrixstructtobasis(m_struct_list):
    matrix_basis_list = []
    for m_struct in m_struct_list:
        matrix_basis_list.append(MatrixStructToBasis(m_struct))
    return matrix_basis_list


def hermitian_csc_matrix(i, j, n):
    if i == j:
        return csc_matrix(([1.0], [[i], [j]]), \
                shape=(n, n), dtype=complex)
    else:
        x = 1.j/np.sqrt(2.)
        return csc_matrix(([x, -x], [[i, j], [j, i]]), \
                shape=(n, n), dtype=complex)


def hermitian_csc_matrix_basis(n, istart=0):
    '''
    Generate hermitian matrix basis set of dimention n.
    '''
    return [hermitian_csc_matrix(i, j, n) for i in range(istart, n) \
            for j in range(i, n)]
