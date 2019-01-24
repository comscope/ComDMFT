from __future__ import print_function

import h5py, numpy, warnings, sys
from scipy.linalg import block_diag
from mpi4py import MPI
from pymatgen.core import Structure
from pymatgen.symmetry.bandstructure import HighSymmKpath
from pymatgen.electronic_structure.bandstructure import BandStructureSymmLine
from pymatgen.electronic_structure.plotter import BSPlotter
from pymatgen.electronic_structure.core import Spin
from builtins import zip,range
import itertools as it
from pyglib.iface.wanniertb import w90


def get_csh2sab():
    '''get transformation from complex spherical harmonics basis to
    g-risb symmetry adapted basis.
    '''
    csh2sab_list = []
    with h5py.File("GPARAM.h5", "r") as f:
        imap_list = f["/IMAP_IMP"][()]
        for i,imap in enumerate(imap_list):
            if i == imap-1:
                csh2sab_list.append(\
                        f["/IMPURITY_{}/DB_TO_SAB".format(i+1)][()].T)
            else:
                csh2sab_list.append(csh2sab_list[imap-1])
    return csh2sab_list


def get_wan2sab():
    '''get transformation from wannier basis to cygutz symmetry adapted basis.
    '''
    csh2sab_list = get_csh2sab()
    with h5py.File("ginit.h5", "r") as f:
        wan2csh = f["/u_wan2csh"][()]
        spin_orb = f["/usrqa/spin_orbit_coup"][()]
    if 'n' in spin_orb.lower():
        iso = 1
    else:
        iso = 2
    csh2sab = []
    for u1 in csh2sab_list:
        n1 = u1.shape[0]/(3-iso)
        csh2sab.append(u1[:n1, ::3-iso])
    csh2sab = block_diag(*csh2sab)
    wan2sab = wan2csh.copy()
    n1 = csh2sab.shape[0]
    wan2sab[:,:n1] = wan2csh[:,:n1].dot(csh2sab)
    return wan2sab


def get_gloc_in_wannier_basis():
    '''get the gutzwiller local matrices, including r, lambda, nr, and nphy
    from grisb calculation.
    '''
    # read from cygutz output
    with h5py.File("GLOG.h5", "r") as f:
        # renormalized local one-body part of the quasiparticle hamiltonian
        lambda_list = f["/BND_LAMBDA"][()].swapaxes(1,2)
        # quasiparticle renormalization matrix
        r_list = f["/BND_R"][()].swapaxes(1,2)
        # local density part to be subtracted in density calculations
        nr_list = f["/BND_NRL"][()].swapaxes(1,2)
        # physical density matrix to be added in density calculations
        nphy_list = f["/BND_NPHY"][()].swapaxes(1,2)

    # get transformation from wannier basis to cygutz symmetry adapted basis.
    wan2sab = get_wan2sab()
    # convert to wannier basis in each spin block.
    lam2 = []
    r2 = []
    nr2 = []
    nphy2 = []
    n2 = wan2sab.shape[0]
    for rmat,lam,nr,nphy in zip(r_list, lambda_list, nr_list, nphy_list):
        n1 = rmat.shape[0]
        if n2 > n1:
            rmat = block_diag(rmat, numpy.eye(n2-n1, dtype=numpy.complex))
            zeromat = numpy.zeros((n2-n1, n2-n1), numpy.complex)
            lam = block_diag(lam, zeromat)
            nr = block_diag(nr, zeromat)
            nphy = block_diag(nphy, zeromat)
        r2.append(wan2sab.dot(rmat).dot(wan2sab.T.conj()))
        lam2.append(wan2sab.dot(lam).dot(wan2sab.T.conj()))
        nr2.append(wan2sab.conj().dot(nr).dot(wan2sab.T))
        nphy2.append(wan2sab.conj().dot(nphy).dot(wan2sab.T))
    return numpy.asarray(r2), numpy.asarray(lam2), \
            numpy.asarray(nr2), numpy.asarray(nphy2)


def get_h1e_in_wannier_basis():
    '''get the h1e-matrix.
    '''
    with h5py.File("GPARAM.h5", "r") as f:
        num_imp = f["/num_imp"][0]
        ispin = f["/ispin"][0]
    # read wannier to csh-basis transformation.
    with h5py.File("ginit.h5", "r") as f:
        wan2csh = f["/u_wan2csh"][()]
    n2 = wan2csh.shape[0]
    h1e_list = []
    with h5py.File("GPARAMBANDS.h5", "r") as f:
        for isp in range(ispin):
            h1e = []
            for i in range(num_imp):
                h1e.append(f["/IMPURITY_{}/H1E_SPIN{}".format(\
                        i+1, isp+1)][()].T)
            h1e = block_diag(*h1e)
            n1 = h1e.shape[0]
            if n2 > n1:
                h1e = block_diag(h1e, numpy.zeros((n2-n1,n2-n1), \
                        dtype=numpy.complex))
            h1e_list.append(wan2csh.dot(h1e).dot(wan2csh.T.conj()))
    return h1e_list


def get_structure():
    with h5py.File("ginit.h5", "r") as f:
        lattice = f["/struct/cell"][()]
        species = f["/struct/symbols"][()]
        coords = f["/struct/scaled_positions"][()]
    return Structure(lattice=lattice, species=species, coords=coords)


def get_bands(kpoints, gmodel=None, wfwannier_list=None,
        bnd_es_in=None, mode="tb"):
    if mode == "risb":
        r_mat, lam_mat, _, _ = get_gloc_in_wannier_basis()
        h1_mat = get_h1e_in_wannier_basis()
        ispin = r_mat.shape[0]
        with h5py.File("GLOG.h5", "r") as f:
            efermi = f["/e_fermi"][0]
    else:
        ispin = 1
        efermi = 0.
    bnd_es = []
    bnd_vs = []
    for isp in range(ispin):
        if mode == "risb":
            ispp = min(isp, len(h1_mat))
        bnd_es.append([])
        bnd_vs.append([])
        for ik, kpt in enumerate(kpoints):
            if gmodel is not None:
                hmat = gmodel._gen_ham(kpt,isp)
            else:
                hmat = wfwannier_list[isp][ik].T.conj().dot(\
                        numpy.diag(bnd_es_in[isp][ik])).dot(\
                        wfwannier_list[isp][ik])
            if mode == "risb":
                hmat -= h1_mat[ispp]
                hmat = r_mat[isp].dot(hmat).dot(r_mat[isp].T.conj())
                hmat += lam_mat[isp]
            evals, evecs = numpy.linalg.eigh(hmat)
            evals -= efermi
            bnd_es[isp].append(evals)
            bnd_vs[isp].append(evecs)
    return numpy.asarray(bnd_es), numpy.asarray(bnd_vs)


def get_symkpath(atol=1.e-6):
    struct = get_structure()
    kpath = HighSymmKpath(struct)
    # check warning and perform transformation if needed.
    if not numpy.allclose(kpath._structure.lattice.matrix,
            kpath._prim.lattice.matrix, atol=atol):
        warnings.warn("Input structure does not match expected standard "
                "primitive! Try k-path transformation.")
        ktrans = kpath._prim.lattice.reciprocal_lattice.matrix.dot(\
                numpy.linalg.inv(kpath._structure.lattice.\
                reciprocal_lattice.matrix))
        for kname in kpath.kpath["kpoints"]:
            kpath.kpath["kpoints"][kname] = \
                    kpath.kpath["kpoints"][kname].dot(ktrans)
    return kpath


def get_gmodel(wpath="../wannier", wprefix="wannier"):
    wannier90 = w90(wpath, wprefix)
    gmodel = wannier90.model()
    return gmodel


def mpiget_bndev(k_list, gmodel=None, wfwannier_list=None, \
        bnd_es_in=None, mode="tb"):
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    ncpu = comm.Get_size()
    nktot = len(k_list)
    nk_cpu = nktot//ncpu
    if nk_cpu*ncpu < nktot:
        nk_cpu += 1
    iksta = rank*nk_cpu
    ikend = min((rank+1)*nk_cpu, nktot)
    kvec_loc = k_list[iksta: ikend]
    if wfwannier_list is not None:
        wfwannier_list = wfwannier_list[:, iksta:ikend, :, :]
        bnd_es_in = bnd_es_in[:, iksta:ikend, :]
    bnd_es, bnd_vs = get_bands(kvec_loc, gmodel=gmodel, \
            wfwannier_list=wfwannier_list, bnd_es_in=bnd_es_in, mode=mode)
    if rank != 0:
        comm.Send(bnd_es, dest=0, tag=rank)
    else:
        for icpu in range(1, ncpu):
            nk_loc = min(nk_cpu, nktot-icpu*nk_cpu)
            _bnd_es = numpy.zeros((bnd_es.shape[0], nk_loc, bnd_es.shape[2]), \
                    dtype=numpy.float)
            comm.Recv(_bnd_es, source=icpu, tag=icpu)
            bnd_es = numpy.concatenate((bnd_es, _bnd_es), axis=1)
        assert bnd_es.shape[1] == nktot, "error in merging bnd_es!"
    return bnd_es, bnd_vs


def get_wannier_den_matrix_risb(bnd_vs, ferwes, wk, nktot):
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    ncpu = comm.Get_size()
    nk_cpu = nktot//ncpu
    if nk_cpu*ncpu < nktot:
        nk_cpu += 1
    r_mat, _, nr_mat, nphy_mat = get_gloc_in_wannier_basis()
    with h5py.File("GPARAMBANDS.h5", "r") as f:
        ispin_dft = f["/ispin"][0]
        iso = f["/iso"][0]
    ispin_risb = r_mat.shape[0]
    # spin factor
    f_ispin = 3-max(iso, ispin_risb)
    wan_den = []
    # total number of electrons to be compared.
    sum_elec1 = numpy.sum(ferwes)
    sum_elec2 = 0.
    for isp in range(ispin_risb):
        if isp <= ispin_dft:
            wan_den.append([])
        else:
            wan_den = numpy.asarray(wan_den)
        for ik, bndvk1, ferwek1, wk1 in zip(it.count(),
                bnd_vs[isp], ferwes[isp], wk):
            # notice the convention a bit different from cygutz.
            # <a|psi>f<psi|b>
            afb = bndvk1.dot(numpy.diag(ferwek1/wk1/f_ispin)).dot(\
                    bndvk1.T.conj())
            # R^\dagger_{A,a} * <a|psi>f<psi|b> * R_{b,B}
            rdafbr = r_mat[isp].T.conj().dot(afb).dot(r_mat[isp])
            dmk = rdafbr
            # \rho_{A,B} = R^\dagger_{A,a} * <a|psi>f<psi|b> * R_{b,B}
            #            +(n_phys.^{A,B} - n_{sub.}^{A,B})
            dmk += (nphy_mat[isp]-nr_mat[isp]).T
            sum_elec2 += dmk.trace()*wk1*f_ispin
            if isp <= ispin_dft:
                wan_den[-1].append(dmk)
            else:
                wan_den[-1][ik] += dmk
                wan_den[-1][ik] *= 0.5
    sum_elec_all1 = comm.reduce(sum_elec1)
    sum_elec_all2 = comm.reduce(sum_elec2)
    if rank == 0:
        elec_diff = sum_elec_all2-sum_elec_all1
        if numpy.abs(elec_diff) > 1.e-4:
            warnings.warn("sum_ferwt = {} vs sum_kswt = {}!". \
                    format(sum_elec1, sum_elec2))
    wan_den = numpy.asarray(wan_den)
    # merge wan_den to master node
    if rank != 0:
        comm.Send(wan_den, dest=0, tag=rank)
    else:
        for icpu in range(1, ncpu):
            nk_loc = min(nk_cpu, nktot-icpu*nk_cpu)
            _wan_den = numpy.zeros((wan_den.shape[0], nk_loc, \
                    wan_den.shape[2], wan_den.shape[3]), dtype=numpy.complex)
            comm.Recv(_wan_den, source=icpu, tag=icpu)
            wan_den = numpy.concatenate((wan_den, _wan_den), axis=1)
        assert wan_den.shape[1] == nktot, "error in merging wan_den!"
    return wan_den


def get_bands_symkpath(efermi=0., mode="tb"):
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    if rank == 0:
        gmodel = get_gmodel()
        kpath = get_symkpath()
        nk = (len(kpath.kpath["kpoints"])-1)*12
        k_vec, k_dist, k_node = gmodel.k_path(kpath.kpath["kpoints"]. \
                values(), nk)
    else:
        gmodel = kpath = k_vec = k_dist = k_node = None
    gmodel = comm.bcast(gmodel, root=0)
    kpath = comm.bcast(kpath, root=0)
    k_vec = comm.bcast(k_vec, root=0)
    k_dist = comm.bcast(k_dist, root=0)
    k_node = comm.bcast(k_node, root=0)
    bnd_es, bnd_vs = mpiget_bndev(k_vec, gmodel=gmodel, mode=mode)
    # prepare the args for pymatgen bs class.
    if rank == 0:
        eigenvals = {}
        eigenvals[Spin.up] = bnd_es[0].T
        if len(bnd_es) == 2:
            eigenvals[Spin.down] = bnd_es[1].T
        bs = BandStructureSymmLine(k_vec, eigenvals, \
                kpath._structure.lattice.reciprocal_lattice, \
                efermi, kpath.kpath["kpoints"])
    else:
        bs = None
    return bs


def plot_bandstructure():
    if "-h" in sys.argv:
        print("usage: complot_bands.py [-g] [-f fname] [-el emin] [-eh emax]")
        sys.exit()

    if "-g" in sys.argv:
        mode = "risb"
    else:
        mode = "tb"
    bs = get_bands_symkpath(mode=mode)
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    if rank == 0:
        bsplot = BSPlotter(bs)
        if "-f" in sys.argv:
            fname = sys.argv[sys.argv.index("-f")+1]
            if ".pdf" not in fname:
                fname += ".pdf"
        else:
            fname = "bndstr.pdf"
        if "-el" in sys.argv:
            emin = float(sys.argv[sys.argv.index("-el")+1])
        else:
            emin = numpy.min(bs.bands.values())
        if "-eh" in sys.argv:
            emax = float(sys.argv[sys.argv.index("-eh")+1])
        else:
            emax = numpy.max(bs.bands.values())
        bsplot.save_plot(fname, img_format="pdf", ylim=(emin, emax), \
                zero_to_efermi=False)
        # better align by yourself, setting zero_to_efermi=False.



if __name__ == "__main__":
    plot_bandstructure()
