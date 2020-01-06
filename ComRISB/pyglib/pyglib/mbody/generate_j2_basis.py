from __future__ import print_function

# generate the j_square basis for each valence block f-shell.
from pyglib.mbody.local_operator_factory import get_nij_op, get_am_op_from_nij
from pyglib.symm.angular_momentum_1p import get_J_vector
from pyglib.symm.unitary import comp_sph_harm_to_relativistic_harm
import numpy, h5py


def h5generate_j2_basis(l="f"):
    '''generate j2 eigen-states with relativistic harmonics as the
    single-particle basis.
    '''
    norbs = {"d":10, "f":14}
    l_list = [(norbs[l]/2-1)/2]
    jvec = get_J_vector(l_list)
    ucsh2rh = comp_sph_harm_to_relativistic_harm(norbs[l])
    jvec = [ucsh2rh.T.conj().dot(jcomp).dot(ucsh2rh) for jcomp in jvec]
    with h5py.File("j2evec.h5", "a") as f:
        group = "/{}".format(l)
        if group in f:
            del f[group]
        for val in range(1,norbs[l]):
            nij = get_nij_op(l=l, ival=val)
            jop_list = get_am_op_from_nij(jvec=jvec, nij=nij,
                    op_list=['Jz', 'J2'])
            j2op = jop_list["J2"].todense()
            jzop = jop_list["Jz"].todense()
            w,v = numpy.linalg.eigh(j2op)
            j_list = map(lambda x: abs(round(numpy.sqrt(x+0.25)-0.5, 1)), w)
            j_unique = numpy.unique(j_list)
            for j in j_unique:
                istr = j_list.index(j)
                iend = len(j_list) - j_list[::-1].index(j)
                v_j = v[:,istr:iend]
                jzop_j = v_j.T.conj().dot(jzop).dot(v_j)
                _w,_v = numpy.linalg.eigh(jzop_j)
                v_j = v_j.dot(_v)
                f["{}/val_{}/j_{}/v".format(group,val,j)] = v_j


def h5calc_dominant_config(l="f", fj2basis="j2evec.h5", \
        frho="EMBED_HAMIL_ANALYSIS_1.h5", fresult="mbbasis.h5", \
        prefix="phy", rcut=1.e-4):
    norbs = {"d":10, "f":14}
    j_ranges = {0:numpy.arange(15), 1:numpy.arange(0.5,15.5)}
    j_list = []
    with h5py.File(fj2basis, "r") as fj:
        with h5py.File(frho, "r") as fin:
            with h5py.File(fresult, "w") as fout:
                for val in range(norbs[l]+1):
                    evals = []
                    evecs = None
                    path_data = "/valence_block_{}/RHO".format(val)
                    if path_data in fin:
                        # fortran to c convention
                        rho = fin[path_data][()].T
                        evals = []
                        if rho.shape[0] == 1:
                            if rho[0,0] > rcut:
                                evals.append(rho[0,0])
                                evecs = [1.0+0.j]
                                j_list.append(0.0)
                        else:
                            for j in j_ranges[val%2]:
                                path_data = "/{}/val_{}/j_{:.1f}/v". \
                                        format(l,val,j)
                                if path_data in fj:
                                    u = fj[path_data][()]
                                    rho_j = -u.T.conj().dot(rho).dot(u)
                                    w,v = numpy.linalg.eigh(rho_j)
                                    w *= -1.
                                    for iw,w1 in enumerate(w):
                                        if w1 < rcut:
                                            if iw == 0 or \
                                                    w[iw-1]-w1 > 1.e-6:
                                                break
                                    if iw < len(w)-1:
                                        iw -= 1
                                    if iw >= 0:
                                        w = w[:iw+1]
                                        v = v[:,:iw+1]
                                        evals.extend(w)
                                        if evecs is None:
                                            evecs = u.dot(v)
                                        else:
                                            evecs = numpy.concatenate(\
                                                    (evecs, u.dot(v)), \
                                                    axis=1)
                                        if (iw+1) % int(round((2*j+1))) != 0:
                                            print("val = {} j = {}".format(\
                                                    val,j))
                                            print(w)
                                            raise ValueError(\
                                                    (" error: j = {} while "+\
                                                    "nevec = {}!").format(\
                                                    j,iw+1))
                                        j_list.extend([j for w1 in w])
                        if evecs is not None:
                            path = "/valence_block_{}".format(val)
                            fout[path+"/{}_rho_evec".format(prefix)] = evecs
                            fout[path+"/{}_rho_eval".format(prefix)] = evals
                if len(j_list) > 0:
                    fout["/{}_j_list".format(prefix)] = j_list



if __name__ == "__main__":
    # h5generate_j2_basis(l="f")
    h5calc_dominant_config()
