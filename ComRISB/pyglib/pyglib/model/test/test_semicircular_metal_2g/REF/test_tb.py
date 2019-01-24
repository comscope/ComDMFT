import numpy as np
import pyglib.gutz.ginput as ginput
import pyglib.model.tbASE as tb
import pyglib.model.special as special


def get_hk1(e, tiny, norb, mu):
    _diag = [e+0.j+mu]
    for i in range(1,norb):
        _diag.append(e*tiny)
    return np.diag(_diag)


def gutz_model_setup(u=0.0, nmesh=5000, norb=1, tiny=0.0, mu=0):
    '''Set up calculations for semi-circular DOS with two-ghost orbitals.
    '''
    # semi-circular energy mesh
    sc = special.semicirular()
    e_list = sc.get_e_list_of_uniform_wt(nmesh=nmesh)

    a = tb.AtomsTB("N", [(0, 0, 0)], cell=(1, 1, 1))
    a.set_orbitals_spindeg(orbitals=[tuple(["s" for i in range(norb)])])

    # set up hk_list
    hk_list = [get_hk1(e, tiny, norb, mu) for e in e_list]
    aTB = tb.TB(a, hk_list=hk_list)

    # dummy k-point list
    kps = e_list
    num_k = len(kps)
    kps_wt = 1.0 / num_k * np.ones((num_k))
    if aTB.Atoms.spindeg:
        kps_wt *= 2
    num_e = 0.0
    norb2 = norb*2

    # GPARAMBANDS.h5

    h1e_list = [np.zeros((norb2, norb2), dtype=np.complex)]
    ginput.save_gparambands(kps_wt, num_e, norb, \
            ensemble=1, h1e_list=h1e_list)

    # GPARAM.h5
    # self-energy structure
    sigma = np.zeros((norb2,norb2), dtype=int)
    sigma[0::2, 0::2] = np.arange(1,norb**2+1).reshape( \
            (norb,norb))
    sigma[1::2, 1::2] = sigma[0::2, 0::2]

    # Coulomb matrix
    v2e = np.zeros((norb2,norb2,norb2,norb2), dtype=np.complex)
    v2e[0,0,0,0] = v2e[0,0,1,1] = v2e[1,1,0,0] = v2e[1,1,1,1] = u
    vdc2_list = np.array([np.zeros((norb2))])
    vdc2_list[0, 0:2] = -u/2

    ginput.save_gparam(na2_list=[norb2],
            iembeddiag=-3, imix=0, sigma_list=[sigma],
            v2e_list=[v2e], nval_top_list=[norb2],
            vdc2_list=vdc2_list)

    # BAREHAM_0.h5
    aTB.save_bareham(kps)


if __name__=='__main__':
    gutz_model_setup(u=2.50, nmesh=5000, norb=3, tiny=0.0)
