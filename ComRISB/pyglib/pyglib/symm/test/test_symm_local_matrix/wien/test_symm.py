from __future__ import print_function
import h5py, numpy
from pyglib.symm.angular_momentum_1p import get_L_vector_CH
from pyglib.symm.atom_symm import get_representation, get_product_table
from pyglib.math.matrix_util import sym_array


# d complex spherical Harmonics
L = get_L_vector_CH(2)

# get Lie parameters and one-particle dm to test
with h5py.File('data_complex_spherical_harmonics.h5', 'r') as f:
    nks = f['/NKS'][()].T
    lie_params = f['/lie_even_params'][()]

# take one-spin component
nks = nks[::2, ::2].T

print('nks.real in comp_sph_harm basis:')
for row in nks:
    print(('{:10.4f}'*len(row)).format(*row.real))

print('nks.imag in comp_sph_harm basis:')
for row in nks:
    print(('{:10.4f}'*len(row)).format(*row.imag))

# rotations
R_list = get_representation(L, lie_params)

# check product table
ptable = get_product_table(R_list)
print('ptable\n',ptable)

for i,rot in enumerate(R_list):
    err = nks.dot(rot) -  rot.dot(nks)
    print('{} maxerr of commutation = {}'.format(i,
            numpy.max(numpy.abs(err))))


nks = numpy.loadtxt('dm_ldau.txt').T

print('nks.real in comp_sph_harm basis:')
for row in nks:
    print(('{:10.4f}'*len(row)).format(*row.real))

print('nks.imag in comp_sph_harm basis:')
for row in nks:
    print(('{:10.4f}'*len(row)).format(*row.imag))

for i,rot in enumerate(R_list):
    err = nks.dot(rot) -  rot.dot(nks)
    print('{} maxerr of commutation = {}'.format(i,
            numpy.max(numpy.abs(err))))

