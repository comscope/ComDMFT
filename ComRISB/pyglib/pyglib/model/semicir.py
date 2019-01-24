# Author: Yongxin Yao

'''One band with semi-circular dos.
'''

import numpy as np
import pyglib.gutz.ginput as ginput
import pyglib.model.tbASE as tb
import pyglib.model.special as special


def gutz_model_setup(u=0.0, nmesh=5000, mu=0):
    '''Set up Gutzwiller calculations for 1-band model with semi-circular DOS.

    Parameters:

    * u: real number
      Hubbard U.
    * nmesh: interger number
      number of energy mesh
    * mu: real number
      chemical potential

    Result:

    Create all the necessary input file of ``GPARAMBANDS.h5``, ``GPARAM.h5``,
    and ``BAREHAM_0.h5``, for *CyGutz* calculation.
    '''
    # get semi-circular class, predifined in pyglib.model.special
    sc = special.semicircular()

    # get semi-circular energy mesh
    e_list = sc.get_e_list_of_uniform_wt(nmesh=nmesh)

    # get atomic structure with 1 atom per unit cell.
    a = tb.AtomsTB("N", [(0, 0, 0)], cell=(1, 1, 1))

    # specify single s-orbital and orbital dimensions.
    a.set_orbitals_spindeg(orbitals=[('s')])
    norb = 1
    norb2 = norb*2

    # get a list of one-body Hamilonian for the enrgy mesh.
    hk_list = [np.array([[e+0.j]]) for e in e_list]

    # get the tight-binding model
    aTB = tb.TB(a, hk_list=hk_list)

    # here the energy mesh is the same as k-points.
    kps = e_list
    num_k = len(kps)

    # set uniform k-point weight.
    kps_wt = 1.0 / num_k * np.ones((num_k))

    # Include the spin degeneracy to k-point weight.
    if aTB.Atoms.spindeg:
        kps_wt *= 2

    # We will work with grand canonical model, num_e will not be used.
    num_e = 0.0

    # list of one-body part of the local Hamiltonians.
    # here the s-orbital is located at zero.
    h1e_list = [np.zeros((norb2, norb2), dtype=np.complex)]

    # generate GPARAMBANDS.h5
    # ensemble is set to 1 for grand canonical system.
    ginput.save_gparambands(kps_wt, num_e, norb, \
            ensemble=1, h1e_list=h1e_list, delta=1.e-8)

    # set self-energy structure
    sigma = np.zeros((norb2,norb2), dtype=int)
    sigma[0::2, 0::2] = np.arange(1,norb**2+1).reshape( \
            (norb,norb))
    sigma[1::2, 1::2] = sigma[0::2, 0::2]

    # set Coulomb matrix
    v2e = np.zeros((norb2,norb2,norb2,norb2), dtype=np.complex)
    v2e[0,0,0,0] = v2e[0,0,1,1] = v2e[1,1,0,0] = v2e[1,1,1,1] = u

    # set the potential shift such that the system has particle-hole
    # symmetry with recpect to zero.
    # also include chemical potential here.
    vdc2_list = np.array([np.zeros((norb2))])
    vdc2_list[0, 0:2] = -u/2 + mu

    # generate GPARAM.h5 file.
    ginput.save_gparam(na2_list=[norb2],
            iembeddiag=-1, imix=0, sigma_list=[sigma],
            v2e_list=[v2e], nval_top_list=[norb2],
            vdc2_list=vdc2_list, max_iter=500)

    # generate BAREHAM_0.h5 file.
    aTB.save_bareham(kps)

