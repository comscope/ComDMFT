import h5py
from mpi4py import MPI
from pyglib.iface.ifwannier import get_wannier_dat
from pyglib.estructure.fermi import get_fermi_level
from pyglib.estructure.gwannier import get_gmodel, mpiget_bndev


kpts, wfwannier_list, bnd_es = get_wannier_dat(path="../wannier")
with h5py.File("GPARAMBANDS.h5", "r") as f:
    num_elec = f['/nelectron'][0]
    ismear = f["/ismear"][0]
    delta = f["/delta"][0]

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
if rank == 0:
    print("num of electrons = {}".format(num_elec))

gmodel = get_gmodel()
kpts = gmodel.k_uniform_mesh((15, 15, 15))
nk = len(kpts)
wklist = [1./nk for k in kpts]
bnd_es, bnd_vs = mpiget_bndev(kpts, gmodel, mode="risb")

if rank == 0:
    print(bnd_es[0][0])
    efermi = get_fermi_level(bnd_es, wklist, num_elec,
            delta=delta, ismear=ismear)
    print("lda fermi level = {}".format(efermi))
