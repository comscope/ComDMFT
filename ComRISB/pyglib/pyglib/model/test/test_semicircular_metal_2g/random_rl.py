import h5py
import numpy as np


def gh5_produce_random_rl(norb, l, z, h):

    norb2 = norb*2

    # R matrix
    r_mat = np.zeros((norb, norb), dtype=complex)
    r_mat[0,0] = z
    r_mat[1,0] = h
    r_mat[2,0] = h
    r_mat2 = np.zeros((norb2, norb2), dtype=complex)
    r_mat2[0::2, 0::2] = r_mat
    r_mat2[1::2, 1::2] = r_mat2[0::2, 0::2]

    # \lambda matrix
    l_mat = np.zeros((norb, norb), dtype=complex)
    l_mat[1,1] = l
    l_mat[2,2] = -l
    l_mat2 = np.zeros((norb2, norb2), dtype=complex)
    l_mat2[0::2, 0::2] = l_mat
    l_mat2[1::2, 1::2] = l_mat2[0::2, 0::2]

    with h5py.File('WH_RL_INIT.h5', 'w') as f:
        f['/IMPURITY_1/LAMBDA'] = l_mat2.T
        f['/IMPURITY_1/R'] = r_mat2.T


if __name__=='__main__':
    gh5_produce_random_rl(norb=3, l=1.1239, z=0.407591572, h=0.6459494)
