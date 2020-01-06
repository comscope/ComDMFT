from __future__ import print_function

import h5py
import numpy as np
from scipy.sparse import csr_matrix
from pyglib.io.h5io import write_csr_matrix
from pyglib.mbody.local_operator_factory import h5get_local_operators
from pyglib.symm.atom_symm import get_representation


def h5get_hs_rotation_representation_N(imp, f_inp="init_ga_info.h5",
        f_log="glog.h5", valences=None):
    '''
    Calculate the Hilbert space rotation representation of each valence block
    of impurity imp in the hilbert_space_rotations.h5 file.
    '''

    f = h5py.File(f_inp, 'r')
    # 0-based in python
    Lie_Jeven_params = f["/impurity_" + str(imp - 1) + "/Lie_Jeven"][...]
    Lie_Jodd_params = f["/impurity_" + str(imp - 1) + "/Lie_Jodd"][...]
    f.close()

    f = h5py.File(f_log, 'r')
    op_dict = h5get_local_operators(imp, f, ["Jx", "Jy", "Jz"])

    rho_nrow = f["/Impurity_" + str(imp) + "/RHO.nrow"][0]
    # one-based to zero-based
    ID_N_block = f["/Impurity_" + str(imp) + "/Sec_ID"][...] - 1
    VAL_N_block = f["/Impurity_" + str(imp) + "/Sec_VAL"][...]
    f.close()
    VAL_N_block = [int(val + 0.5) for val in VAL_N_block]
    if valences is None:
        valences = VAL_N_block

    f = h5py.File("hilbert_space_rotations.h5", 'a')
    base = op_dict["Jx"].shape[0] - rho_nrow

    for i, val in enumerate(VAL_N_block):
        if ID_N_block[i + 1] - ID_N_block[i] <= 0:
            continue
        if val not in valences:
            continue
        Jx = op_dict["Jx"][base + ID_N_block[i] : base + ID_N_block[i + 1],
                           base + ID_N_block[i] : base + ID_N_block[i + 1]]
        Jy = op_dict["Jy"][base + ID_N_block[i] : base + ID_N_block[i + 1],
                           base + ID_N_block[i] : base + ID_N_block[i + 1]]
        Jz = op_dict["Jz"][base + ID_N_block[i] : base + ID_N_block[i + 1],
                           base + ID_N_block[i] : base + ID_N_block[i + 1]]
        if abs(np.mod(val, 2)) > 1.e-4:
            Rpr = get_representation([Jx, Jy, Jz], Lie_Jodd_params)
        else:
            Rpr = get_representation([Jx, Jy, Jz], Lie_Jeven_params)

        f["Impurity_{}/val_block={}/dim_rotations".format(imp, val)] = \
                len(Rpr)

        for j, _R in enumerate(Rpr):
            _R = csr_matrix(_R)
            write_csr_matrix(f,
                    "Impurity_{}/val_block={}/rotation_{}".format(imp, val, j),
                    _R)
    f.close()


if __name__ == "__main__":
    h5get_hs_rotation_representation_N(1, valences=[2])
