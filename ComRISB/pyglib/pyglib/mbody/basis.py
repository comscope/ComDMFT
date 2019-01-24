from __future__ import print_function
try:
    from builtins import range
except:
    pass

import h5py, numpy


def get_full_fock_states_by_int(norb):
    '''generate all the fock states represented by integers.
    '''
    fock_list = [[] for i in range(norb+1)]
    for i in range(2**norb):
        n1 = bin(i).count('1')
        fock_list[n1].append(i)
    return fock_list


def get_double_fock_states_at_half_filling(fock_list):
    '''generate the impurity-bath fock states at half-filling.
    '''
    norb = len(fock_list)-1
    base = 2**norb
    nstates = sum([len(fs)**2 for fs in fock_list])
    fock2_array = numpy.empty(nstates, dtype=numpy.int32)
    iadd = 0
    for i, fs in enumerate(fock_list):
        for fs1 in fs:
            for fs2 in fock_list[norb-i]:
                fock2_array[iadd] = fs1 + fs2*base
                iadd += 1
    return fock2_array


def h5add_fock_states(norb, f_h5):
    '''add the fock state basis.
    '''
    fock_list = get_full_fock_states_by_int(norb)
    fock2_array = get_double_fock_states_at_half_filling(fock_list)
    f_h5["/fock_basis/dim"] = [len(fock2_array)]
    f_h5["/fock_basis/vals"] = fock2_array


if __name__ == "__main__":
    with h5py.File("fock.h5", "w") as f:
        h5add_fock_states(14, f)
