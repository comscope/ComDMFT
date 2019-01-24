from __future__ import print_function

import numpy as np
import h5py
from scipy.linalg import block_diag
from pyglib.io.fio import file_exists


def get_orbital_transformation(l, nspin):
    '''give the unitary transformation from wien2k basis to the standard
    complex spherical Harmonics with Condon-shortley phase.
    '''

    u_trans = np.identity((2*l+1), dtype=complex)

    # Test for the weird phase (identified in NiO AFM case.)
    if l == 2:
        u_trans[1,1] = -1.j
        u_trans[3,3] =  1.j

    if nspin == 2:
        u_trans = block_diag(u_trans, u_trans)
    return u_trans


def h4set_indmfl(emin=-10., emax=10., projtype=-2):
    '''create case.indmfl for the interface given file ginit.h5.
    '''
    finit = h5py.File('ginit.h5', 'r')
    if '/struct/case' not in finit:
        raise ValueError('path /struct/case does not exist in ginit.h5!')

    case = finit['/struct/case'][()].split('.')[0]

    if file_exists(case+'.indmfl'):
        print(case+'.indmfl exists! \n' +
                'Please rm it if you intend to generate a new one.')
        return

    findm = open(case+'.indmfl', 'w')
    findm.write(('{:6.2f} {:6.2f} {:3d}  # hybrid. emin/max w.r.t. FS,' + \
            ' projector type\n').format(emin, emax, projtype))

    atom_symbols = finit['/struct/symbols'][()]
    unique_corr_atom_symbols = finit['/usrqa/unique_corr_symbol_list'][()]. \
            tolist()
    n_corr_atoms = 0
    for symbol in atom_symbols:
        if symbol in unique_corr_atom_symbols:
            n_corr_atoms += 1
    findm.write('{:2d}                 # number of correlated atoms\n'. \
            format(n_corr_atoms))

    if 'y' in finit['/usrqa/spin_orbit_coup'][()]:
        nspin = 2
    else:
        nspin = 1
    unique_df_list = finit['/usrqa/unique_df_list'][()]

    from pyglib.symm.angular_momentum_1p import get_l_list_from_string
    cix = 0
    norbsmax = 0
    l_super_list = []
    for i, symbol in enumerate(atom_symbols):
        if symbol not in unique_corr_atom_symbols:
            continue
        df = unique_df_list[unique_corr_atom_symbols.index(symbol)]
        l_list = get_l_list_from_string(df)
        findm.write('{:2d} {:2d}              # iatom, num_L\n'. \
                format(i+1, len(l_list)))
        l_super_list.extend(l_list)
        for l in l_list:
            cix += 1
            findm.write(' {:2d} {:2d} {:2d}          # L, qsplit, cix\n'.\
                    format(l, 3, cix))
            norbsmax = max(norbsmax, (2*l+1)*nspin)

    findm.write('----- unitary transformation for correlated orbitals -----\n')
    findm.write(('{:2d} {:2d}              # num indep. kcix blocks,' + \
            ' max dimension\n').format(cix, norbsmax))
    for i, l in enumerate(l_super_list):
        norbs = (2*l+1)*nspin
        findm.write('{:2d} {:2d}              # cix, dimension\n'. \
                format(i+1, norbs))
        findm.write('----- unitary transformation matrix -----\n')
        u_trans = get_orbital_transformation(l, nspin)
        for u1 in u_trans:
            for u12 in u1:
                findm.write('{:20.16f}{:20.16f} '.format(u12.real, u12.imag))
            findm.write('\n')


if __name__=='__main__':
    h4set_indmfl()
