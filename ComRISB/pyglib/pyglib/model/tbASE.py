# Author: Yongxin Yao, Xiaoyu Deng

'''
A tight-binding model is defined by a set of hoppoing matrix H(R).
H(R) is a matrix with dimension of orbitals in the unit cell.
Dimension: all cells are three dimension in ASE,
however one can mimick two dimensional model
by setting hopping along the third axis zero.
'''

import numpy
import h5py
from ase import Atoms
from ase.dft import kpoints


class AtomsTB(Atoms):

    """
    ASE ATOMS object with orbitals.

    Attributes:

      - For attributes of ASE ATOMS object,
        please check ASE documentation.
      - orbitals: orbital names/index,
        list with each item is a list of orbital names
        on each atom, for example, [["s"],["s","p"]]
      - spin: name of spins, by default, spin=["up"],
        only one spin component is considered.
        Can be also ["up","dn"].
      - nspinorbitals: number of orbitals with possible
        spin degeneracy
    """

    def set_orbitals_spindeg(self, orbitals=None,
            spindeg=True, spinorbit=False):
        """
        Set the orbitals in a ASE Atoms object.

        Parameter:
        ----------
          orbitals: list
            orbitals names, list with each item
            is a list of orbital names on each atom,
            for example, [("s"),("s","px")]
          spindeg: bool
            spin degeneracy, if TRUE spin=["up"];
            if False, spin=["up","dn"]
          spinorbit: bool
            True, orbitals are in fact spin-orbitals,
            if this is true, spin is set to ["so"],
            False, orbitals are without spinorbit coupling.

        Example::

        >>> # set a atom with two bands, "s" and "p",
        >>> # and spin degenaracy
        >>> a=AtomsTB("N",[(0,0,0)],cell=(1,1,1))
        >>> a.set_orbitals_spindeg(orbitals=[("s","px")],
        >>>         spindeg=True)
        """

        if orbitals is None:
            # one s orbital per atoms by default.
            self.orbitals = [("s",) for pos in self.positions]
        else:
            self.orbitals = orbitals[:]

        self.spindeg = spindeg
        self.spinorbit = spinorbit
        if spinorbit:
            self.nso = 2
        else:
            self.nso = 1

        self.nspinorbitals = 0
        for i in self.orbitals:
            self.nspinorbitals += len(i)
        self.nspinorbitals *= self.nso


class TB(object):

    """
    A tight-binding mode is defined with an ASE object Atoms
    and Hopping.
    """

    def __init__(self, AtomsTB, Hr=None, hk_list=None):
        """
        Init TB object

        Parameters
        ----------
        AtomsTB: ASE object
          define the unit cell, the atoms, with orbitals set.
        Hr: dict
          a dict of (R, hop), with R is tranlational vector
          and hop is hopping matrix between unit cells.
          len(Hr[R])=self.AtomsTB.nspinorbitals.
        """

        self.Atoms = AtomsTB
        self.Hr = Hr
        self.hk_list = hk_list


    @staticmethod
    def gallery(name="SimpleSquareLattice"):
        """
        A gallery of predefined lattice that frequentely used.

        Parameters
        ----------
        name: model name
          currently "SimpleSquareLattice", "Chain_nn" implemented.

        Returns
        -------
        aTB: TB object
        """

        if name == "SimpleSquareLattice":
            a = AtomsTB("N", [(0, 0, 0)], cell=(1, 1, 1))
            a.set_orbitals_spindeg()
            aTB = TB(a)
            aTB.set_hop([((0, 1, 0), 0, 0, -1),
                         ((1, 0, 0), 0, 0, -1),
                         ((0, -1, 0), 0, 0, -1),
                         ((-1, 0, 0), 0, 0, -1)])
            return aTB

        if name == "Chain_nn":
            a = AtomsTB("N", [(0, 0, 0)], cell=(1, 1, 1))
            a.set_orbitals_spindeg()
            aTB = TB(a)
            aTB.set_hop([((1, 0, 0), 0, 0, 1),
                         ((-1, 0, 0), 0, 0, 1)])
            return aTB

    def set_hop(self, hoppings=None):
        """
        Set hopping mannually. In addition to a Hr matrix.
        Hopping is a tupe with the form (R,iorb,jorb,t_hop).
        hoppings could be a list of hoppings.

        Parameters
        ----------
        hoppings: list
          a list of hopping defined as tuple(R,iorb,jorb,t_hop).
          For example, [((0,0,0),0,0,1)], or a tuple(R,hr),
          where hr is a matrix of hoppings

        Example::

        >>>a=AtomsTB("N",[(0,0,0)],cell=(1,1,1))
        >>>a.set_orbitals_spindeg()
        >>>a.set_hop([((0,1,0),0,0,-1),
        >>>           ((1,0,0),0,0,-1)])
        """
        # if only set only one term, change it to list for consistency.
        if type(hoppings) == type((1,)):
            hoplist = [hoppings, ]
        else:
            hoplist = hoppings
        if self.Hr is None:
            self.Hr = {}
        for ihop in hoplist:             #
            if len(ihop) > 2:
                R, iorb, jorb, t_hop = ihop
                if R not in self.Hr:  # the matrix for R is not set yet
                    self.Hr[R] = numpy.zeros((
                            self.Atoms.nspinorbitals,
                            self.Atoms.nspinorbitals), dtype=numpy.complex)
                self.Hr[R][iorb, jorb] = t_hop
            else:
                R, hr = ihop
                self.Hr[R] = hr


    def Hk(self, ikp=0, kpt=numpy.zeros((3))):
        """
        Construct Hamiltonian for a given k point, i.e.,
        fourier transformation of Hr.

        Parameters
        ----------
        kpt: array
          size (3). k point is required to be in unit of reciprocal
          periodic lattice vector.

        Returns
        -------
        Hk: ndarray
          dimension (norbitals,norbitals), with each element
          the Hamiltonian for the corresponding kpoint.
        kcart: ndarray
          k points in cartesian coordinates.
        """

        if self.Hr is not None:
            reciprocalvector = self.Atoms.get_reciprocal_cell()
            self.hk = numpy.zeros((self.Atoms.nspinorbitals, \
                    self.Atoms.nspinorbitals), dtype=numpy.complex)
            ikcart = numpy.dot(kpt, reciprocalvector)
            for ir in self.Hr:
                R = numpy.array(ir)
                Rcart = numpy.dot(R, numpy.array(self.Atoms.cell))
                expk = numpy.exp(-1j * ikcart.dot(Rcart) * 2.0 * numpy.pi)
                self.hk += expk * self.Hr[ir]
        elif self.hk_list is not None:
            self.hk = self.hk_list[ikp]
            ikcart = None
        else:
            raise ValueError('Neither Hr nor hk_list are defined!')

        return self.hk, ikcart

    def eigens(self, ikp=0, kpt=numpy.zeros((3))):
        """
        Caculate eigenvalues of a given list of kpoints.
        It calls Hk and diagonalize the hamiltonian.

        Parameters
        ----------
        kpt: numpy.ndarry
          kpoint, dimension (3), in unit of reciprocal lattice vectors.

        Returns
        -------
        eks: ndarray
          eigenvalues, taken to be real.
        Uk: ndarray
          wavefunctions.
        """

        hk, _ = self.Hk(ikp, kpt)
        if hk.shape[0] == 1:
            ek, Uk = [hk[0,0].real], numpy.array([[1.0+0.j]])
        else:
            ek, Uk = numpy.linalg.eigh(hk)
        return ek, Uk

    def get_bandstructure(self, kps, saveto=None, with_weights=False):
        """
        Get the band structure for given kpath.

        Parameters
        ----------
        kps: list
          kpoints, kpath object given by ase.dft.kpoints.
        saveto: string
          name of the file to save data.
        with_weights:

        """

        eks = []
        if with_weights:
            Uks = []
        for i, kpt in enumerate(kps[0]):
            ek, Uk = self.eigens(i, kpt)
            eks.append(ek)
            if with_weights:
                Uks.append(Uk)
        if saveto is None:
            filename = "band.dat"
        else:
            filename = saveto
        with open(filename, "w") as f:
            for iorb in range(self.Atoms.nspinorbitals):
                for ik in range(len(kps[0])):
                    f.write("%f  %f " % (kps[1][ik], eks[ik][iorb]))
                    if with_weights:
                        for ilay in range(self.Atoms.nspinorbitals):
                            weight = Uks[ik, ilay, iorb] * \
                                Uks[ik, ilay, iorb].conj()
                            f.write("%f  " % (weight.real))
                    f.write("\n")
                f.write("\n")

    def get_dos(self,  kps_size=None, saveto=None, dos_mesh=None, eta=1e-2):
        """
        Get the density of states. Note: symmetry is not used.

        Parameters
        ----------
        dos_mesh: ndarray
          the mesh of dos. default: numpy.linspace(Emax,Emin,500),
          (Emax, Emin)= (max,min) among all the eigenvalues
        kps_size: tuple
          the size of Monkhorst_pack mesh. default: (8,8,8)
        saveto: string
          name of the file to save data.
        eta: float
          broadening factor. current only Gaussian braodening available,
          default: 1e-2

        Returns
        -------
        dos: density of states

        Examples::

        >>># dos of 2D square lattice
        >>>aTB=TB.gallery()
        >>>aTB.get_dos(kps_size=(400,400,1))
        """

        if kps_size is None:
            kps_size = (8, 8, 8)
        kps = kpoints.monkhorst_pack(kps_size)
        eks = [self.eigens(i, kpt)[0] for i, kpt in enumerate(kps)]

        if dos_mesh is None:
            dos_mesh = numpy.linspace(eks.min(), eks.max(), 500)
        if saveto is None:
            filename = "dos.dat"
        else:
            filename = saveto
        dos = numpy.zeros(len(dos_mesh), dtype=numpy.float)
        for ek in numpy.nditer(eks):
            dos += numpy.exp(-(ek - dos_mesh)**2 / 2.0 / eta**2)
        dos *= 1.0 / numpy.sqrt(2 * numpy.pi) / eta
        dos /= len(kps)

        with open(filename, "w") as f:
            for idx in xrange(len(dos_mesh)):
                f.write("%f  %f \n" % (dos_mesh[idx], dos[idx]))

    def add_spindegeneracy(self):
        """
        Add spin degeneracy to a non-spin TB hamiltonian.
        Note simply increase the number of orbitals and increase the size
        of Hr to 2x2 block diagonal form. up and down spin
        in block-diagonal form.

        Returns
        -------
        TB: AtomsTB
          A new AtomsTB object with spin degeneracy.
        """
        atoms = self.Atoms.copy()
        assert not self.Atoms.spindeg, \
                " Error: spin degeneray is already considered!"
        assert not self.Atoms.spinorbit, \
                " Error: Cann't add spin degeneracy to spin-orbit orbitals!"
        atoms.set_orbitals_spindeg(self.Atoms.orbitals, spindeg=True)
        norb = atoms.nspinorbitals
        if self.Hr is not None:
            Hr = {}
            for iR in self.Hr:
                Hr[iR] = numpy.zeros((norb, norb), \
                        dtype=type(self.Hr[iR][0, 0]))
                Hr[iR][0:norb / 2, 0:norb / 2] = self.Hr[iR][:, :]
                Hr[iR][norb / 2:, norb / 2:] = self.Hr[iR][:, :]
        else:
            Hr = None
        return TB(atoms, Hr)


    def save_bareham(self, kpts, ngroup=1):
        nk_per_group = len(kpts)/ngroup + 1
        ik_start_list = [i*nk_per_group for i in range(ngroup)]
        for igroup, ik_start in enumerate(ik_start_list):
            with h5py.File('BAREHAM_'+str(igroup)+'.h5', 'w') as f:
                ik_end = min(ik_start+nk_per_group, len(kpts))
                for ik in range(ik_start, ik_end):
                    ek, Uk = self.eigens(ik, kpts[ik])
                    # Upper case: Fortran convention
                    f['/IKP_'+str(ik+1)+'/ek0'] = ek
                    f['/IKP_'+str(ik+1)+'/ISYM_1/HK0'] = self.hk.T
                    f['/IKP_'+str(ik+1)+'/T_PSIK0_TO_HK0_BASIS'] = Uk.T



if __name__ == "__main__":
    # The following is a simple test for the above codes.

    # AtomsTB object
    # set up an AtomsTB object. one band without spin degeneracy.
    a = AtomsTB("N", [(0, 0, 0)], cell=(1, 1, 1))
    a.set_orbitals_spindeg()

    # sqare lattice
    # set up a TB object and the hoppings. This corresponds to a 2D squared
    # lattice with nearest-neighbor hopping. This is the default one in
    # TB.gallery()
    aTB = TB(a)
    aTB.set_hop([((0, 1, 0), 0, 0, -1),
                 ((1, 0, 0), 0, 0, -1),
                 ((0, -1, 0), 0, 0, -1),
                 ((-1, 0, 0), 0, 0, -1)])

    # bands and dos
    # set special k points for bands
    kG = [0, 0, 0]
    kX = [0.5, 0, 0]
    kM = [0.5, 0.5, 0]
    # set up a ase.dft.kpoints kpath object
    kps = kpoints.get_bandpath([kG, kX, kM, kG], a.cell)
    # get band structure of a square lattice
    aTB.get_bandstructure(kps, saveto="pcell_band.dat")
