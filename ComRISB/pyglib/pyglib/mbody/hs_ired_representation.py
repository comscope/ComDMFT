from __future__ import print_function


import h5py
from pyglib.io.h5io import get_csr_matrix
from pyglib.symm.atom_symm import get_atom_U_sigma_chi


def get_U_chi_sigma(imp, valences):
    '''
    Get the U_list, chi_list and sigma_list in the Hilbert space.
    '''
    f = h5py.File("hilbert_space_rotations.h5", 'r')
    list_U = []
    list_chi = []
    list_sigma = []
    for val in valences:
        Rpr = []
        for j in range(f["Impurity_{}/val_block={}/dim_rotations".format( \
                imp, val)][...]):
            Rpr.append(get_csr_matrix(f, \
                    "Impurity_{}/val_block={}/rotation_{}".format( \
                    imp, val, j)).todense())
        U, sigma, chi = get_atom_U_sigma_chi(Rpr)
        list_U.append(U)
        list_chi.append(chi)
        list_sigma.append(sigma)
    f.close()
    return list_U, list_chi, list_sigma


if __name__ == "__main__":
    list_U, list_chi, list_sigma = get_U_chi_sigma(1, [2])
    print(list_sigma[0])
    for chi in list_chi[0]:
        print(''.join(" {:<4.1f}".format(x) for x in chi))
    print(list_U)
