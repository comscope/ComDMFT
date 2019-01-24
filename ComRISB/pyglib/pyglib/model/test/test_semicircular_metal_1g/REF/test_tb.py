import numpy as np
import pyglib.gutz.ginput as ginput
import pyglib.model.tbASE as tb
import pyglib.model.special as special


def gutz_model_setup(u=0.0):
    '''Set up calculations for semi-circular DOS with two-ghost orbitals.
    '''
    # semi-circular energy mesh
    sc = special.semicirular()
    e_list = sc.get_e_list_of_uniform_wt()

    a = tb.AtomsTB("N", [(0, 0, 0)], cell=(1, 1, 1))
    a.set_orbitals_spindeg(orbitals=[("p")])
    # set up hk_list
    hk_list = [np.diag([0., e - u/2, 0.j]) for e in e_list]
    aTB = tb.TB(a, hk_list=hk_list)

    # dummy k-point list
    kps = e_list
    num_k = len(kps)
    kps_wt = 1.0 / num_k * np.ones((num_k))
    if aTB.Atoms.spindeg:
        kps_wt *= 2
    num_e = 0.0
    num_band_max = 3

    # GPARAMBANDS.h5
    h1e_list = [np.array([[0,0,0],[0, -u/2,0],[0,0,0]], \
            dtype=np.complex)]
    ginput.save_gparambands(kps_wt, num_e, num_band_max, \
            ensemble=1, h1e_list=h1e_list)

    # GPARAM.h5
    # self-energy structure
    sigma = np.zeros((6,6), dtype=int)
    sigma[0::2, 0::2] = np.arange(1,10).reshape((3,3))
    sigma[1::2, 1::2] = sigma[0::2, 0::2]

    # Coulomb matrix
    v2e = np.zeros((6, 6, 6, 6), dtype=np.complex)
    v2e[2, 2, 2, 2] = v2e[2, 2, 3, 3] = v2e[3, 3, 2, 2] = v2e[3, 3, 3, 3] = u

    ginput.save_gparam(na2_list=[6], iembeddiag=-1, imix=0, sigma_list=[sigma],
            v2e_list=[v2e], nval_top_list=[6])

    # BAREHAM_0.h5
    aTB.save_bareham(kps)


if __name__=='__main__':
    gutz_model_setup(u=1.0)
