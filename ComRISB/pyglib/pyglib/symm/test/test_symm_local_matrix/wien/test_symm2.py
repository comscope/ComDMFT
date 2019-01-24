from __future__ import print_function
import h5py, numpy
from pyglib.symm.angular_momentum_1p import get_L_vector_CH
from pyglib.symm.atom_symm import get_representation, get_product_table
from pyglib.math.matrix_util import sym_array


# d complex spherical Harmonics
L = get_L_vector_CH(2)

# get Lie parameters and one-particle dm to test
with h5py.File('data_complex_spherical_harmonics.h5', 'r') as f:
    lie_params = f['/lie_even_params'][()]

# rotations
R_list = get_representation(L, lie_params)

# check product table
ptable = get_product_table(R_list)
print('ptable\n',ptable)


with h5py.File('FeSb2_local_matrix.h5', 'r') as f:
    nks_sang = f['/NKS_PRE_SYM_Sang'][()].T[::2,::2]
    nks_wien = f['/NKS_PRE_SYM_wien'][()].T[::2,::2]

for i,rot in enumerate(R_list):
    err = nks_sang.dot(rot) -  rot.dot(nks_sang)
    print('{} maxerr of commutation = {}'.format(i,
            numpy.max(numpy.abs(err))))

for i,rot in enumerate(R_list):
    err = nks_wien.dot(rot) -  rot.dot(nks_wien)
    print('{} maxerr of commutation = {}'.format(i,
            numpy.max(numpy.abs(err))))

