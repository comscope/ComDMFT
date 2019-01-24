import os, sys, subprocess, h5py, numpy, shlex
from pyglib.gutz import embedh
from pyglib.mbody.basis import h5add_fock_states


def create_ged_inp(imp, h1e, lambdac, daalpha, v2e, istep=0, rtol=1.e-10,
        nev=1, maxiter=1000):
    '''create/update the input file for the parallel ed solver.
    '''
    with h5py.File("GEDINP_{}.h5".format(imp), "a") as f:
        if v2e is not None:
            ind_list = numpy.where(numpy.abs(v2e) > rtol)
            val_list = v2e[ind_list]
            # in fortran convention, one-based
            f["/v2e/nnz"] = [len(val_list)]
            ind_list = numpy.array(ind_list)+1
            f["/v2e/INDS"] = ind_list.T
            # absorb the 1/2 factor.
            f["/v2e/vals"] = val_list/2.

        v1e = embedh.get_whole_h1e(h1e, lambdac, daalpha)
        ind_list = numpy.where(numpy.abs(v1e) > rtol)
        val_list = v1e[ind_list]
        if "/v1e" in f:
            del f["/v1e"]
        norb = v1e.shape[0]
        f["/v1e/norb"] = [norb]
        f["/v1e/nnz"] = [len(val_list)]
        ind_list = numpy.array(ind_list)+1
        f["/v1e/INDS"] = ind_list.T
        f["/v1e/vals"] = val_list

        if "/ed" not in f:
            f["/ed/nev"] = [nev]
            f["/ed/maxiter"] = [maxiter]
            f["/ed/rtol"] = [rtol]
            h5add_fock_states(norb/2, f)


def driver_ed(imp=1, istep=0, mpiexec="mpirun -np 2 "):
    '''dmrg solver based on ITensor.
    '''
    if "-i" in sys.argv:
        imp = int(sys.argv[sys.argv.index("-i")+1])
    if os.path.isfile("GEDINP_{}.h5".format(imp)):
        lv2e = False
    else:
        lv2e = True
    h1e, lambdac, daalpha, v2e = embedh.h5gen_embedh_spin_updn(
            imp=imp, lv2e=lv2e)
    create_ged_inp(imp, h1e, lambdac, daalpha, v2e, istep=istep)

    if os.path.isfile("mpi_ed.txt"):
        mpiexec = open("mpi_ed.txt", "r").readline()

    mpiexec = shlex.split(mpiexec)
    cmd = mpiexec+[os.environ["WIEN_GUTZ_ROOT2"]+"/exe_ed", "-i", str(imp)]

    print(" running {}".format(" ".join(cmd)))
    subprocess.call(cmd)

    # get the calculation results in cmp_sph_harm basis with
    # faster spin index.
    with h5py.File("GEDOUT_{}.h5".format(imp)) as f:
        dm = f["/DM"][()].T
        emol = f["/emol"][0]

    # transform it to the symmetry-adpated basis
    embedh.h5wrt_dm_sab_cmplx(dm, emol, imp=imp)



if __name__ == "__main__":
    driver_ed()
