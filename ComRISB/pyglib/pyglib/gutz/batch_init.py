from __future__ import print_function

import h5py
from pyglib.gutz.init import initialize


def batch_initialize(cell, scaled_positions, symbols, case=None,
        crystal_field='y', dist_cut=-1, full_orbital_polarization='y',
        idx_equivalent_atoms=None, iembeddiag=-1, ldc=1, lnewton=0,
        spin_orbit_coup='n', spin_polarization='n', u_matrix_type=1,
        unique_corr_symbol_list=None, unique_df_list=None,
        unique_j_list_ev=None, unique_nf_list=None,
        updn_full_list=None,
        unique_u_list_ev=None, unit='eV'):
    '''initialize *CyGutz* calculation by directly providing values
    of the list of arguments.

    Parameters:

    * cell: 3x3 matrix
        Unit cell vectors.
    * scaled_positions: ist of scaled-positions
        Atomic positions given in units of the unit cell.
    * symbols: list of str
        a list of symbols. E.g., ['H', 'H', 'O']
    * case: str
        case name for Wien2k.
    * crystal_field: str
        'y' for considering crystal field, and 'n' for not
    * dist_cut: real number
        cut-off distance for extracting finite fragment centered at an atom
        from the infinite lattice. If default (-1), the code would determine
        it automatically.
    * full_orbital_polarization: str
        y' for breaking all the possible orbital symmetry, and 'n' for not.
    * idx_equivalent_atoms: list of integers
        Indices for equivalent atoms. E.g., 0 0 0 1 1' means 1-3 and 4-5
        are two inequivalent atoms.
    * iembeddiag: integer
        flag for the method to solve embedding Hamiltonian.

        * -3: Valence truncation ED for S=0 (spin-singlet)
        * -1: Valence truncation ED
        * 10: Hartree-Fock (LDA+U)

    * ldc: integer
        flag for Coulomb interaction double counting.

        * 12: Recommended for LDA+G-RISB calculations.
            FLL double counting.
            (updating Vdc at each charge iteration,
            initial n0 to be provided.)
        *  1: FLL double counting potential.
            (n_0 self-consistently determined.)
        *  0: No double counting (useful for models).

    * lnewton: integer
        flag for the method to solve G-RISB equations.

        *  0: Modified Powell hybrid method (HYDRD1).
        * -1: Broyden method.

    * spin_orbit_coup: str
        'y' for considering spin-orbit coupling, and 'n' for not.
    * spin_polarization: str
        'y' for spin-polarization, and 'n' for not
    * u_matrix_type: integer
         flag for the parametrize Coulomb U-matrix.

         * 1: Slater-Condo parametrization.
         * 2: Kanamori parametrization (useful for models).
         * 0: Manual input.

    * unique_corr_symbol_list: list of str
        list of unique correlated atom symbols
    * unique_df_list: list of str
        list of s/p/d/f orbitals
        for the list of unique correlated atom symbols
    * unique_j_list_ev: list of real
        list of Hund's parameter J
        for the list of unique correlated atom symbols
    * unique_u_list_ev: list of real
        list of Hubbard U
        for the list of unique correlated atom symbols
    * unique_nf_list: list of real
        list of initial impurity orbital occupations
        for the list of unique correlated atom symbols
    * unit: str
        label for the unit system: 'eV' or 'rydberg'.

    Results:

    finish the initializztion for *CyGutz* accordingly and create
    the ``GPARAM.h5`` file.
    '''

    with h5py.File('ginit.h5', 'w') as f:
        if case is not None:
            f['/struct/case'] = case
        f['/struct/cell'] = cell
        f['/struct/scaled_positions'] = scaled_positions
        f['/struct/symbols'] = symbols
        f['/usrqa/crystal_field'] = crystal_field
        f['/usrqa/dist_cut'] = dist_cut
        f['/usrqa/full_orbital_polarization'] = full_orbital_polarization
        if idx_equivalent_atoms is None:
            idx_equivalent_atoms = range(len(symbols))
        f['/usrqa/idx_equivalent_atoms'] = idx_equivalent_atoms
        f['/usrqa/iembeddiag'] = iembeddiag
        f['/usrqa/ldc'] = ldc
        f['/usrqa/lnewton'] = lnewton
        f['/usrqa/spin_orbit_coup'] = spin_orbit_coup
        f['/usrqa/spin_polarization'] = spin_polarization
        if 'y' in spin_polarization:
            if updn_full_list is None:
                updn_full_list = [1 for s in symbols]
            f['/usrqa/updn_full_list'] = updn_full_list
        f['/usrqa/u_matrix_type'] = u_matrix_type
        if unique_corr_symbol_list is None:
            unique_corr_symbol_list = list(set(symbols))
        f['/usrqa/unique_corr_symbol_list'] = unique_corr_symbol_list
        if unique_df_list is None:
            unique_df_list = ['s' for x in unique_corr_symbol_list]
        f['/usrqa/unique_df_list'] = unique_df_list
        if unique_j_list_ev is None:
            unique_j_list_ev = [0. for x in unique_corr_symbol_list]
        f['/usrqa/unique_j_list_ev'] = unique_j_list_ev
        if unique_u_list_ev is None:
            unique_u_list_ev = [0. for x in unique_corr_symbol_list]
        f['/usrqa/unique_u_list_ev'] = unique_u_list_ev
        if unique_nf_list is None:
            unique_nf_list = [1. for x in unique_corr_symbol_list]
        f['/usrqa/unique_nf_list'] = unique_nf_list
        f['/usrqa/unit'] = unit

    initialize()
