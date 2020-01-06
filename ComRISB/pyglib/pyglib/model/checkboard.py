import numpy as np
import pickle
from ase.dft import kpoints
import pyglib.gutz.ginput as ginput
import pyglib.model.tbASE as tb

def gutz_model_setup(u=0.0, spindeg=True, num_e=2., iembeddiag=-1):
    '''Set up Gutzwiller calculations for 2d body-centered square lattice.

    Parameters:

    * u: real number
      Hubbard U.
    * spindeg: boolean number
      whether to keep spin degeneracy or not.
    * num_e: real number
      number of electron per unit cell
    * iembeddiag: integer
      flag for method to solve the embedding Hamiltonian.

      * -3: valence truncation ED with S=0 (spin-singlet) constraint;
      * -1: valence truncation ED;
      * 10: Hartree-Fock.

    Result:

    Create all the necessary input file of ``GPARAMBANDS.h5``, ``GPARAM.h5``,
    and ``BAREHAM_0.h5``, for *CyGutz* calculation.
    '''

    # two "H" atoms in the square unit cell, one at the corner
    # and one at the center.
    symbols=['H', 'H']
    scaled_positions=[(0, 0, 0), (0.5, 0.5, 0)]
    cell = np.identity(3)
    a = tb.AtomsTB(symbols=symbols, scaled_positions=scaled_positions,
            cell=cell)

    # set spin degeneracy accordingly.
    a.set_orbitals_spindeg(spindeg=spindeg)

    # create a tight-binding model class given the AtomsTB.
    aTB = tb.TB(a)

    # set real space (nearest neighbour) hopping elements.
    t = -1.0
    aTB.set_hop([
            (( 0, 0,0),0,1, t),
            ((-1, 0,0),0,1, t),
            (( 0,-1,0),0,1, t),
            ((-1,-1,0),0,1, t),
            (( 0, 0,0),1,0, t),
            (( 1, 0,0),1,0, t),
            (( 0, 1,0),1,0, t),
            (( 1, 1,0),1,0, t),
            ])

    # set 2d k-mesh
    kps_size = (50, 50, 1)
    kps = kpoints.monkhorst_pack(kps_size)

    # set uniform k-point weight
    num_k = len(kps)
    kps_wt = 1.0 / num_k * np.ones((num_k))
    if aTB.Atoms.spindeg:
        kps_wt *= 2

    # se maximal number of bands (here we have two bands.)
    num_band_max = 2

    # set list of one-body part of the local Hamiltonian (trivial here.)
    h1e_list = [np.array([[0,]], dtype=np.complex) for symbol in symbols]

    # create ``GPARAMBANDS.h5`` file
    ginput.save_gparambands(kps_wt, num_e, num_band_max, h1e_list=h1e_list)

    # cerate ``BAREHAM_0.h5`` file.
    aTB.save_bareham(kps)

    # for initializing CyGutz calculation
    from pyglib.gutz.batch_init import batch_initialize

    if spindeg:
        spin_polarization = 'n'
    else:
        spin_polarization = 'y'

    # the default settings are good for this model.
    batch_initialize(cell=cell, scaled_positions=scaled_positions,
            symbols=symbols,idx_equivalent_atoms=[0,1],
            unique_u_list_ev=[u], iembeddiag=iembeddiag,
            spin_polarization=spin_polarization,
            updn_full_list=[1,-1])

    # save the aTB
    with open('aTB.pckl', 'wb') as f:
        pickle.dump([a, aTB], f)


def get_bard_bands():
    '''get bare band structure.
    '''
    # load aTB
    with open('aTB.pckl', 'rb') as f:
        a, aTB = pickle.load(f)

    # k-point path
    kG = [0.0, 0.0, 0]
    kX = [0.5, 0.0, 0]
    kM = [0.5, 0.5, 0]

    # set up a ase.dft.kpoints kpath object
    kps = kpoints.get_bandpath([kG, kX, kM, kG], a.cell)

    # get band structure of a square lattice
    aTB.get_bandstructure(kps, saveto="bare_bands.dat")


if __name__=='__main__':
    gutz_model_setup()
    get_bard_bands()
