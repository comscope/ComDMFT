import numpy as np
from scipy.linalg import block_diag
from pyglib.math.matrix_util import trans_orbital_fast_to_spin_fast
from pyglib.symm.angular_momentum_1p import get_J_generator, \
        get_JU_relat_sph_harm_cg
import pyglib.symm.atom_symm as atsym


def get_self_energy(l_list, ispin, orbital_pol, cf, iso, rotations=None):
    '''
    Dispatcher.
    '''
    if 'y' in orbital_pol and iso == 2:
        return get_self_energy_op_soc(l_list)
    elif 'y' in cf and iso == 2:
        return get_self_energy_cf_soc(l_list, rotations)
    elif iso == 2:
        return get_self_energy_soc(l_list, ispin)
    elif 'y' in orbital_pol:
        return get_self_energy_op_nosoc(l_list, ispin)
    elif 'y' in cf:
        return get_self_energy_cf_nosoc(l_list, ispin, rotations)
    else:
        return get_self_energy_average(l_list, ispin)


def get_self_energy_op_soc(l_list):
    '''
    Get the self energy structure in the case of fully orbital symmetry
    breaking and with spin-orbit interaction.
    '''
    dim_tot = np.sum(2 * (2 * np.array(l_list) + 1))
    self_energy = np.arange(dim_tot * dim_tot, dtype=int) + 1
    self_energy = self_energy.reshape((dim_tot, dim_tot))
    U = np.identity(dim_tot, dtype=complex)
    return None, U, self_energy


def get_self_energy_op_nosoc(l_list, ispin):
    '''
    Get the self energy structure in the case of fully orbital symmetry
    breaking but with negligible spin-orbit interaction.
    '''
    dim_t = int(np.sum(2 * np.array(l_list) + 1) + 0.5)
    m_half = np.arange(dim_t * dim_t, dtype=int) + 1
    m_half = m_half.reshape((dim_t, dim_t))
    if ispin == 2:
        shift = np.max(m_half)
    else:
        shift = 0
    from scipy.linalg import block_diag
    self_energy = block_diag(m_half, m_half + shift)

    # spin-fast-index convention
    self_energy = trans_orbital_fast_to_spin_fast(self_energy)
    Uhalf = np.identity(dim_t)

    # {{orbs}_up, {orbs}_dn} -> {{up,dn}_{orbs}}
    U = np.zeros((dim_t*2, dim_t*2), dtype=complex)
    U[:dim_t,0::2] = Uhalf
    U[dim_t:,1::2] = Uhalf
    return None, U, self_energy


def get_self_energy_cf_soc(l_list, rotations):
    '''
    Get the self energy structure in the case of crystal field splitting
    and with spin-orbit interaction.
    '''
    j_rel, u_csh2rel = get_JU_relat_sph_harm_cg(l_list)
    J, U, self_energy = atsym.get_atom_Jnew(rotations, j_rel)
    U = u_csh2rel.dot(U)
    return J, U, self_energy


def get_self_energy_cf_nosoc(l_list, ispin, rotations):
    '''
    Get the self energy structure in the case of crystal field splitting
    and without spin-orbit interaction.
    '''
    Jorig = get_J_generator(l_list, iso=1)
    Jhalf, Uhalf, self_energy_half = atsym.get_atom_Jnew(rotations, Jorig)

    if ispin == 2:
        shift = np.max(self_energy_half)
    else:
        shift = 0

    self_energy = block_diag(
        self_energy_half, shift_self_energy(self_energy_half, shift))

    # spin-fast-index convfention
    self_energy = trans_orbital_fast_to_spin_fast(self_energy)

    dim_t = Uhalf.shape[0]
    J = []
    for Jh1 in Jhalf:
        J1 = np.zeros(self_energy.shape, dtype=np.complex)
        J1[0::2,0::2] = Jh1
        J1[1::2,1::2] = Jh1
        J.append(J1)

    U = np.zeros(self_energy.shape, dtype=np.complex)
    U[:dim_t,0::2] = Uhalf
    U[dim_t:,1::2] = Uhalf
    return J, U, self_energy


def get_self_energy_average(l_list, ispin):
    '''
    Get the self energy structure in the case of negligible crystal field
    splitting and without spin-orbit interaction.
    '''
    diag_elem = []
    elem_base = 1
    for l in l_list:
        diag_elem += [elem_base for i in range(2 * l + 1)]
        elem_base += 1
    elem_base -= 1
    if ispin == 2:
        diag_elem += [i + elem_base for i in diag_elem]
    else:
        diag_elem += diag_elem
    self_energy = np.diag(diag_elem)
    self_energy = trans_orbital_fast_to_spin_fast(self_energy)
    U = np.identity(len(self_energy), dtype=complex)
    n_half = U.shape[0]/2
    utrans = np.zeros_like(U)

    # Similarly transform to spin-fast convention.
    utrans[:,::2] = U[:,:n_half]
    utrans[:,1::2] = utrans[:,::2]
    return None, utrans, self_energy


def get_self_energy_soc(l_list, ispin):
    '''
    Get the self energy structure in the case that spin-orbit interaction
    is dominant.
    '''
    diag_elem = []
    elem_base = 1
    for l in l_list:
        if l == 0:
            if ispin == 1:
                elem = [elem_base, elem_base]
            else:
                elem = [elem_base, elem_base + 1]
        else:
            elem = []
            if ispin == 1:
                elem += [elem_base for i in \
                        range(int(2 * (l - 0.49)) + 1)]
                elem_base = max(elem) + 1
                elem += [elem_base for i in \
                        range(int(2 * (l + 0.51)) + 1)]
            else:
                elem += [elem_base + i for i in \
                        range(int(2 * (l - 0.49)) + 1)]
                elem_base = max(elem) + 1
                elem += [elem_base + i for i in \
                        range(int(2 * (l + 0.51)) + 1)]
        diag_elem += elem
        elem_base = max(elem) + 1
    self_energy = np.diag(diag_elem)

    J, U = get_JU_relat_sph_harm_cg(l_list)
    return J, U, self_energy


def shift_self_energy(self_energy, shift):
    '''
    Shift non-zero self-energy elements by shift.
    '''
    res = np.zeros_like(self_energy)
    spots = np.where(self_energy > 0.1)
    res[spots] = self_energy[spots] + shift
    return res



if __name__ == "__main__":
    '''
    A test.
    '''
    l_list = [1, 1]
    ispin = 1
    orbital_pol = 'n'
    cf = 'y'
    iso = 1
    rotations = [[[1, 0, 0], [0, 1, 0], [0, 0, 1]]]
    U, self_energy = get_self_energy(
        l_list, ispin, orbital_pol, cf, iso, rotations)
    print " self_energy = "
    print self_energy

    # H cluster
    symbols = ['H', 'H', 'H', 'H', 'H', 'H']
    positions = [[0., 0., 0.], [2.0, 0., 0.], [0., 2.0, 0.0],
                 [0., -2.0, 0.], [0., 0., 2.0], [0., 0., -2.0]]
    from pyglib.gutz.molecule import xyz_get_rot_list
    rotations = xyz_get_rot_list(symbols, positions, log='screen')
    print "H-dimer rotations:"
    print rotations
    l_list = [0, 0, 1]
    orbital_pol = 'n'
    cf = 'y'
    J, U, self_energy = get_self_energy(
        l_list, ispin, orbital_pol, cf, iso, rotations)
    print " self_energy = "
    print self_energy
