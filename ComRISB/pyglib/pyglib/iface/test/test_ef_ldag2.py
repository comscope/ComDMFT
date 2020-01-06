import h5py, numpy
from mpi4py import MPI
from pyglib.estructure.fermi import get_fermi_level
from pyglib.estructure.gwannier import get_gmodel, mpiget_bndev
from pyglib.iface.ifwannier import get_risb_bndes


with h5py.File("GPARAMBANDS.h5", "r") as f:
    num_elec = f['/nelectron'][0]
    ismear = f["/ismear"][0]
    delta = f["/delta"][0]
    kpts = f["/kpoints"][()]

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
if rank == 0:
    print("num of electrons = {}".format(num_elec))

gmodel = get_gmodel()
#kpts = gmodel.k_uniform_mesh((6, 6, 6))
nk = len(kpts)
wklist = [1./nk for k in kpts]
bnd_es, bnd_vs = mpiget_bndev(kpts, gmodel, mode="risb")
if rank == 0:
    bnd_es2 = get_risb_bndes()
    for ik in range(nk):
        if numpy.max(numpy.abs(bnd_es2[0][ik] - bnd_es[0][ik])) > 1.e-6:
            print(ik)
            print(bnd_es2[0][ik])
            print(bnd_es[0][ik])
    efermi = get_fermi_level(bnd_es, wklist, num_elec,
            delta=delta, ismear=ismear)
    print("lda fermi level = {}".format(efermi))
