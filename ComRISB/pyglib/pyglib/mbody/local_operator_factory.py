'''
Get lcoal operators given the {n_ij} operators.
'''

import h5py, pkg_resources
import itertools as it
import numpy as np
from scipy.sparse import csc_matrix
from pyglib.io.h5io import get_coo_matrix


def get_linear_operators(nij, uij):
    '''
    Calculate \sum_{i,j} nij_{i,j}*uij(i,j).
    '''
    if nij is None:
        return None
    res = None
    for i, n1, u1 in it.izip(it.count(), nij, uij):
        for j, n12, u12 in it.izip(it.count(), n1, u1):
            if n12 is None:
                continue
            if abs(u12)>1.e-12:
                if res is None:
                    res = u12*n12
                else:
                    res += u12*n12
            if i==j:
                continue
            if abs(uij[j,i])>1.e-12:
                if res is None:
                    res = uij[j,i]*n12.getH()
                else:
                    res += uij[j,i]*n12.getH()
    if res is None:
        res = csc_matrix((1,1))
    return res


def get_am_op_from_nij(svec=None, lvec=None, jvec=None,
        nij=None, op_list=None):
    '''
    Get angular momentum operators given the {n_ij} operators and
    coiefficient matrice of spin and orbital momentum operators.
    '''
    if op_list is None:
        return None

    if svec is not None:
        s_z = get_linear_operators(nij, svec[2])
        s_p = get_linear_operators(nij, svec[0]+1.j*svec[1])
    if lvec is not None:
        l_z = get_linear_operators(nij, lvec[2])
        l_p = get_linear_operators(nij, lvec[0]+1.j*lvec[1])
    if svec is not None and lvec is not None:
        j_z, j_p = s_z + l_z, s_p + l_p
    elif jvec is not None:
        j_z = get_linear_operators(nij, jvec[2])
        j_p = get_linear_operators(nij, jvec[0]+1.j*jvec[1])

    res = {}

    if "Sx" in op_list:
        res["Sx"] = (s_p.getH() + s_p) / 2.0
    if "Sy" in op_list:
        res["Sy"] = (s_p - s_p.getH()) / 2.0j
    if "Sz" in op_list:
        res["Sz"] = s_z
    if "Lx" in op_list:
        res["Lx"] = (l_p.getH() + l_p) / 2.0
    if "Ly" in op_list:
        res["Ly"] = (l_p - l_p.getH()) / 2.0j
    if "Lz" in op_list:
        res["Lz"] = l_z
    if "Jx" in op_list:
        res["Jx"] = (j_p.getH() + j_p) / 2.0
    if "Jy" in op_list:
        res["Jy"] = (j_p - j_p.getH()) / 2.0j
    if "Jz" in op_list:
        res["Jz"] = j_z
    if "S2" in op_list:
        res["S2"] = s_p.getH() * s_p + s_z * s_z + s_z
    if "L2" in op_list:
        res["L2"] = l_p.getH() * l_p + l_z * l_z + l_z
    if "J2" in op_list:
        res["J2"] = j_p.getH() * j_p + j_z * j_z + j_z
    for key in res.keys():
        res[key] = res[key].tocsc()
    return res


def get_local_operators(imp, ival, op_list):
    '''
    Get the local operators, like Sx, Sy, Sz, Lx, Ly, Lz, S^2, L^2
    and J^2 operators for impurity imp.
    '''
    # Read in coefficient matrices for S and L.
    with h5py.File('GPARAM.h5', 'r') as f:
        svec = []
        svec.append(f['/IMPURITY_'+str(imp)+'/SX'][...].T)
        svec.append(f['/IMPURITY_'+str(imp)+'/SY'][...].T)
        svec.append(f['/IMPURITY_'+str(imp)+'/SZ'][...].T)
        lvec = []
        lvec.append(f['/IMPURITY_'+str(imp)+'/LX'][...].T)
        lvec.append(f['/IMPURITY_'+str(imp)+'/LY'][...].T)
        lvec.append(f['/IMPURITY_'+str(imp)+'/LZ'][...].T)

    # Read n_ij operators
    norb = svec[0].shape[0]
    nij = []
    with h5py.File('EMBED_HAMIL_ANALYSIS_'+str(imp)+'.h5', 'r') as f:
        for i in range(norb):
            n_i = []
            for j in range(i+1):
                path = '/valence_block_{}/NP_{}_{}'.format(ival,i+1,j+1)
                if path in f:
                    n_i.append(get_coo_matrix(f, path).tocsc())
                else:
                    n_i.append(None)
            nij.append(n_i)

    return get_am_op_from_nij(svec=svec, lvec=lvec, nij=nij, op_list=op_list)


def get_label_list(J, vecs):
    '''
    Get list of expection values of J for vecs.
    '''
    if J is None:
        label = None
    else:
        label = []
        for vec in vecs:
            res = np.vdot(vec, J.dot(vec))
            assert np.abs(res.imag) < 1.e-8, " error in get_label_list!"
            label.append(res.real)
        label = np.array(label)
    return label


def get_nij_op(l='d', ival='1'):
    '''get the matrix representation of n_{i,j}=c_{i}^{\dagger} c_{j}
    in the Hilbert space with valence of ival from preset file 'mbody.h5'.
    '''
    if l == 'd':
        norb = 10
    elif l == 'f':
        norb = 14
    else:
        raise ValueError('unsupported orbital={}!'.format(l))
    fname = pkg_resources.resource_filename('pyglib', './mbody/mbody.h5')
    with h5py.File(fname, 'r') as f:
        nij = []
        for i in range(norb):
            n_i = []
            for j in range(i+1):
                path = '/{}/valence_block_{}/n_{}_{}'.format(
                        l,ival,i,j)
                n_i.append(get_coo_matrix(f, path).tocsc())
            nij.append(n_i)
    return nij


def get_lvec_op(lvec, l='d', ival='1'):
    '''get the l angular momentum operators [lx, ly, lz]
    in the Hilbert space with valence of ival,
    given their representation in single particle space (lvec).
    '''
    nij = get_nij_op(l=l, ival=ival)
    lops = get_am_op_from_nij(lvec=lvec, nij=nij,
            op_list=['Lx', 'Ly', 'Lz'])
    return [lops['Lx'], lops['Ly'], lops['Lz']]


def get_lrot_op(lvec, lie_jeven_params, l='d', ival='1'):
    '''get the 3d rotation matrix representations
    in the Hilbert space with valence of ival,
    given the single particle space representation of the l-angular momentum
    operators and lie parameters.
    '''
    lops = get_lvec_op(lvec, l=l, ival=ival)
    from pyglib.symm.atom_symm import get_representation
    reps = get_representation(lops, lie_jeven_params)
    return reps
