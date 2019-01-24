from __future__ import print_function
'''
Help functions for io.
'''
try:
    from builtins import range
except:
    pass

import numpy as np
import h5py, glob


def read_slice(lines):
    n = len(lines) / 2
    U = []
    for i in range(n):
        V = []
        line_real = lines[i]
        line_real = line_real.split()
        line_imag = lines[i + n]
        line_imag = line_imag.split()
        for j in range(n):
            V.append(float(line_real[j]) + float(line_imag[j]) * 1.j)
        U.append(V)
    U = np.array(U).T  # Fortran convension to C convension
    return U


def read_TRANS(fileName):
    '''
    Read transformation matrix.
    '''
    with open(fileName, 'r') as f:
        lines = f.readlines()
    line = lines[0]
    line = line.split()
    numOrbitals = int(line[0])
    numIons = int(line[1])
    u_list = []
    iBase = 1
    for ion in range(numIons):
        u_list.append(read_slice(lines[iBase: iBase + numOrbitals * 2]))
        iBase += numOrbitals * 2
    return u_list


def read_Hs(fileName):
    '''
    Read Hermitian matrix basis.
    '''
    with open(fileName, 'r') as f:
        lines = f.readlines()
    line = lines[0]
    line = line.split()
    numIons = int(line[0])
    iLine = 1
    HsAll = []
    for ion in range(numIons):
        Hs = []
        line = lines[iLine]
        line = line.split()
        numOrbitals = int(line[0])
        numHs = int(line[1])
        iLine += 1
        for iHs in range(numHs):
            Hs.append(read_slice(lines[iLine: iLine + numOrbitals * 2]))
            iLine += numOrbitals * 2
        HsAll.append(Hs)
    return HsAll


def write_slice(f, U):
    for V in U.T:
        f.write("".join("%20.16f" % (x.real) for x in V) + '\n')
    for V in U.T:
        f.write("".join("%20.16f" % (x.imag) for x in V) + '\n')


def write_TRANS(fileName, u_list):
    '''
    Write transformation matrix.
    '''
    numIons = len(u_list)
    numOrbitals = len(u_list[0])
    f = open(fileName, 'w')
    f.write(" %5i %5i" % (numOrbitals, numIons) + '\n')
    for u in u_list:
        write_slice(f, u)
    f.close()


def write_RLNEF(fileName, RAll, LA1All, NKSAll, EFermi):
    '''
    Write new solution under the current basis transformation.
    '''
    numIons = len(RAll)
    numOrbitals = len(RAll[0])
    with open(fileName, 'w') as f:
        f.write(" %5i %5i" % (numOrbitals, numIons) + '\n')
        for ion in range(numIons):
            write_slice(f, RAll[ion])
            write_slice(f, LA1All[ion])
            write_slice(f, NKSAll[ion])
        f.write(" %20.12f" % (EFermi))


def write_Hs(fileName, HsAll):
    '''
    Write new solution under the current basis transformation.
    '''
    numIons = len(HsAll)
    with open(fileName, 'w') as f:
        f.write(" %5i" % (numIons) + '\n')
        for ion in range(numIons):
            numHs = len(HsAll[ion])
            f.write(" %5i %5i" % (len(HsAll[ion][0]), numHs) + '\n')
            for iHs in range(numHs):
                write_slice(f, HsAll[ion][iHs])


def write_sigma_struct(sigma_list, file_name="WH_SIGMA_STRUCT.INP"):
    '''
    Write WH_SIGMA_STRUCT.INP file.
    '''
    with open(file_name, 'w') as f:
        for ni, sigma in enumerate(sigma_list):
            f.write(" %5i %5i" % (ni + 1, len(sigma)) + '\n')
            for col in np.array(sigma).T:
                f.write(' '.join('%4d' % (x) for x in col) + '\n')


def write_frozen(ne_mott_list, orb_mott_list):
    '''
    Write FROZEN.INP file.
    '''
    with open("FROZEN.INP", 'w') as f:
        for ni, orbs in enumerate(orb_mott_list):
            f.write(" %5i ! NI" % (ni + 1) + '\n')
            f.write(" %5i %5i ! N_ORB, N_ELECTRON_MOTT_LOCALIZED" % (
                    len(orbs), ne_mott_list[ni]) + '\n')
            if len(orbs) > 0:
                f.write(''.join(" %5i " % (orb + 1) for orb in orbs) + '\n')


def write_WH_SL_VEC(s_vec_list, l_vec_list):
    '''
    Write WH_SL_VEC.INP file.
    '''
    s_vec_list = np.swapaxes(s_vec_list, 0, 1)
    l_vec_list = np.swapaxes(l_vec_list, 0, 1)
    with open("WH_SL_VEC.INP", 'w') as f:
        for _s_list in s_vec_list:
            for _s in _s_list:
                write_slice(f, _s)
        for _l_list in l_vec_list:
            for _l in _l_list:
                write_slice(f, _l)


def special_reorder(a_list):
    if len(a_list) != 1 or a_list[0].shape[0] != 6:
        print('Dimension not correct, special reorder cancelled!')
        return a_list

    umat = np.array([
        [0, 1, 0],
        [1, 0, 0],
        [0, 0, 1]], dtype=complex)

    mat2 = np.zeros((6, 6), dtype=complex)
    mat2[0:3, 0:3] = umat.dot(a_list[0].T[::2, ::2]).dot(umat.T.conj()).T
    mat2[3:6, 3:6] = mat2[0:3, 0:3]
    return [mat2]


def h52txt_rln(fh5_name, ftxt_name, specialt=False):
    lma_list = []
    nks_list = []
    r_list = []
    with h5py.File(fh5_name, 'r') as f:
        for imp in range(1000):
            path = '/IMPURITY_{}'.format(imp+1)
            if path in f:
                lma_list.append(f[path+'/LAMBDA'][()].T)
                nks_list.append(f[path+'/NKS'][()].T)
                r_list.append(f[path+'/R'][()].T)
            else:
                break
    if specialt:
        lma_list = special_reorder(lma_list)
        nks_list = special_reorder(nks_list)
        r_list = special_reorder(r_list)

    write_RLNEF(ftxt_name, r_list, lma_list, nks_list, 0.0)


def h52txt_rln_list(h5p='WH_RL_INIT.h5', txtp='WH_RLNEF.INP', specialt=False):
    for fname in glob.glob(h5p+'*'):
        h52txt_rln(fname, txtp+fname[len(h5p):], specialt=specialt)



if __name__=='__main__':
    h52txt_rln_list(specialt=True)
