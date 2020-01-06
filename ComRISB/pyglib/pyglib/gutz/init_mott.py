#!/usr/bin/env python

import h5py
import numpy as np
from pyglib.math.matrix_basis import sigmatomatrixbasis, matrixstructtobasis
from pyglib.gutz.usrqa import get_usr_input

def get_shrink_sigma(sigma):
    '''
    Shrink the structure of sigma if numbers skipped.
    '''
    elem_max = np.max(sigma)
    for elem in range(elem_max):
        if elem == np.max(sigma):
            break
        while elem + 1 not in sigma:
            spots = np.where(sigma > elem)
            sigma[spots] -= 1
    return sigma


def print_sigma(sigma, title):
    '''
    Print the structure of sigma with basis indices.
    '''
    print(title)
    print("index  " + ''.join("%4d " % (i) for i in range(len(sigma))) + '\n')
    for i, row in enumerate(sigma):
        print("%4d   " % (i) + ''.join("%4d " % (j) for j in row))


def init_mott(fin='GPARAM.h5', fout='GMOTT.h5'):
    '''
    Generate the input files for Mott phase calculation.
    '''
    dim_hs_r_list = []
    dim_hs_l_list = []
    f = h5py.File(fin, 'r')
    fp = h5py.File(fout, 'w')
    num_imp = f["/num_imp"][...]
    for imp in range(num_imp):
        prepath = "/IMPURITY_" + str(imp+1)
        sigma = f[prepath + "/SIGMA_STRUCT"][...]
        print "**********  Impurity " + str(imp+1) + "  **********"
        print_sigma(sigma, ' Sigma structure:')
        while True:
            orb_mott = raw_input(
                    ' Please provide the indices of orbitals to'
                    +' be Mott localized \n (e.g., 0 2 ): ')
            yn = raw_input(
                ' You selected [' +
                orb_mott +
                '] \n to be Mott localized, right? (y/n):')
            if 'y' in yn or 'Y' in yn:
                orb_mott = orb_mott.split()
                ind_orb_mott = [int(s) for s in orb_mott]
                while True:
                    ne_mott = raw_input(
                            ' Please provide the total number of'
                            +' Mott localized electrons (per unit cell): ')
                    yn = raw_input(
                        ' Total ' +
                        ne_mott +
                        ' electrons will be Mott localized, right? (y/n):')
                    if 'y' in yn or 'Y' in yn:
                        break
                break
        fp[prepath+'/num_mott_orbitals'] = [len(ind_orb_mott)]
        fp[prepath+'/num_mott_electrons'] = [int(ne_mott)]
        fp[prepath+'/mott_orbital_indices'] = np.asarray(ind_orb_mott)+1

        for ind in ind_orb_mott:
            sigma[ind, :] = 0
        sigma = get_shrink_sigma(sigma)
        print_sigma(sigma, ' R structure:')
        hs = matrixstructtobasis(sigma)
        dim_hs_r_list.append(len(hs))

        fp[prepath+'/SIGMA_STRUCT_R'] = sigma.T
        if len(hs) > 0:
            fp[prepath+'/HS_R'] = np.swapaxes(hs,1,2)

        for ind in ind_orb_mott:
            sigma[:, ind] = 0
        sigma = get_shrink_sigma(sigma)
        print_sigma(sigma, ' Lambda structure:')
        hs = sigmatomatrixbasis(sigma)
        dim_hs_l_list.append(len(hs))
        fp[prepath+'/SIGMA_STRUCT_L'] = sigma.T
        if len(hs) > 0:
            fp[prepath+'/HS_L'] = np.swapaxes(hs,1,2)

    print("\n Please select the method to solve embedding Hamiltonian.\n" +
            " LDIAG = -2: Valence truncation ED for Mott solution.\n" +
            "         -3: Option (-2) with additional S=0 constraints.\n" +
            "         -4: Option (-2) with additional J=0 constraints.\n" +
            "        -31: Option (-2) with S=0 \n" +
            "             and local symmetry constraints.")
    iembeddiag = get_usr_input(" Please select LDIAG: ",
            ['-2','-3','-4','-31'])

    fp['/giembeddiag'] = [int(iembeddiag)]
    fp['/dim_hs_r_imp'] = dim_hs_r_list
    fp['/dim_hs_l_imp'] = dim_hs_l_list
    f.close()
    fp.close()


if __name__ == "__main__":
    init_mott()
