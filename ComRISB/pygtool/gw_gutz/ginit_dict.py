cygutz_init={
'/struct/cell': [\
        [5.8328, 0, 0], \
        [0, 6.5376, 0], \
        [0, 0, 3.1973]],  # [a1, a2, a3]

'/struct/scaled_positions': [\
        [0, 0, 0],
        [0.5, 0.5, 0.5],
        [0.1882, 0.3552, 0],
        [0.8118, 0.6448, 0],
        [0.6882, 0.1448, 0.5],
        [0.3118, 0.8552, 0.5]], # [pos_atom1, pos_atom_2, ...]

'/usrqa/unit': 'eV',
'/usrqa/spin_polarization': 'n',
'/usrqa/full_orbital_polarization': 'y', # do not consider local symmetry
'/usrqa/spin_orbit_coup': 'n',
'/usrqa/crystal_field': 'y', # must since full_orbital_polarization = 'y'

'/usrqa/ldc': 1, # fully localized limit DC
'/usrqa/iembeddiag': -1, # Exact diagonalization
'/usrqa/lnewton': 0, # newton solver
'/usrqa/u_matrix_type': 1, # Slater-condon parametrization

'/struct/symbols': ["Fe", "Fe", "Sb", "Sb", "Sb", "Sb"],
'/usrqa/idx_equivalent_atoms': [0, 0, 2, 2, 2, 2], # Fe1 is the same as Fe0
'/usrqa/unique_corr_symbol_list': ["Fe"],
'/usrqa/unique_df_list': ["d"],
'/usrqa/unique_u_list_ev': [4.8],
'/usrqa/unique_j_list_ev': [0.9]
}


# covert it to hdf5 format.
import h5py

with h5py.File('ginit.h5', 'w') as f:
    for key in cygutz_init:
        f[key] = cygutz_init[key]
