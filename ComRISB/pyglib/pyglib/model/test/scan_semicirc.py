import os, sys, h5py
import numpy as np
from pyglib.model import semicir as bethlatt
import subprocess
from builtins import zip

import matplotlib
matplotlib.rcParams['backend'] = "Qt5Agg"
import matplotlib.pyplot as plt


def generate_data(u_list, mu_list):
    '''run *CyGutz* calculations for a list of U and mu.
    save the results in results.

    Parameters:

    * u_list: real array
      list of Hubbard U parameters.
    * mu_list: real array
      list of chemical potential

    Result:

    save the list of energies (e_list), double occupancy (d_list),
    quasi-particle weight (z_list) and occupation (n_list)
    to ``result.dat`` text file as well as the ``result.h5`` hdf5 file.
    '''
    # get *CyGutz* command. Choose option '-r -1' to skip the claculation of
    # Gutzwiller renormalized electron density.
    root = os.environ['WIEN_GUTZ_ROOT2']
    cmd = [root+'/CyGutz', '-r', '-1']

    # set chemical potential to 0.
    mu=0

    # set Hubbard U=0
    u = 0.

    # generate input files with updated with u=0, mu=0
    bethlatt.gutz_model_setup(u=u, nmesh=5000, mu=mu)

    # total energy list
    e_list = []

    # quasi-particle weight z list
    z_list = []

    # double occupancy d list
    d_list = []

    # occupation list
    n_list = []

    # loop over the list of U.
    for u, mu in zip(u_list, mu_list):

        print(' working on u = {} mu = {}'.format(u, mu))

        # modify the local Coulomb matrix
        with h5py.File('GPARAM.h5', 'a') as f:

            # note the transposition, which is transformation
            # from Fortran convention to c-convention.
            v2e = f['/IMPURITY_1/V2E'][()].T
            vdc2_list = f['VDC2_LIST'][()].T

            # now update the Coulom matrix
            v2e[0,0,0,0] = v2e[0,0,1,1] = v2e[1,1,0,0] = v2e[1,1,1,1] = u
            f['/IMPURITY_1/V2E'][()] = v2e.T

            # update vdc2, which keeps the particle-hole symmetry of the model
            # for finite U.
            vdc2_list[0, 0:2] = -u/2 + mu
            f['VDC2_LIST'][()] = vdc2_list.T

        # perform the *CyGutz* calculation.
        subprocess.call(cmd)

        # get total energy
        with h5py.File('GLOG.h5', 'r') as f:
            nks = f['/IMPURITY_1/NKS_SYM'][()]
            n = nks[0, 0] + nks[1, 1]
            n_list.append(n.real)
            e = f['/etot_model'][0] + u/2*n -mu*n
            e_list.append(e.real)

        # get Z = R^\dagger R
        with h5py.File('WH_RL_OUT.h5', 'r') as f:
            r = f['/IMPURITY_1/R'][0,0]

            # simple here since it is just a scalar
            z = r*r.conj()
            z_list.append(z.real)

        # to get double occupancy (of impurity 1), <n_up n_dn>_G,
        # we run analysis code *exe_spci_analysis*
        subprocess.call([root+'/exe_spci_analysis', '1'])

        # double occupancy is simply the local many-body density matrix element
        # in the valence=2 block.
        with h5py.File('EMBED_HAMIL_ANALYSIS_1.h5', 'r') as f:
            if '/valence_block_2/RHO' in f:
                d = f['/valence_block_2/RHO'][0, 0]
            else:
                d = 0.
            d_list.append(d.real)

    with open('result.dat', 'w') as f:
        for u, mu, e, z, d, n in zip(u_list, mu_list,
                e_list, z_list, d_list, n_list):
            f.write('{} {} {} {} {} {}\n'.format(u, mu, e, z, d, n))

    with h5py.File('result.h5', 'w') as f:
        f['/mu_list'] = mu_list
        f['/u_list'] = u_list
        f['/e_list'] = e_list
        f['/z_list'] = z_list
        f['/d_list'] = d_list
        f['/n_list'] = n_list


def scan_u(mu=0.0):
    '''run *CyGutz* calculations for a list of U for a given mu.

    Parameters:

    * mu: real number
      the fixed chemical potential.

    Result:

    it will generate results for a u_list of np.arange(0.0, 5.1, 0.2)
    at fixed mu.
    '''
    if os.path.isfile('result.h5'):
        return

    # set range of Hubbard U.
    u_list = np.arange(0.0, 5.1, 0.2)
    mu_list = [mu for u in u_list]
    generate_data(u_list, mu_list)


def scan_mu(u=5.0):
    '''run *CyGutz* calculations for a list of mu for a given u.

    Parameters:

    * u: real number
      the fixed Hubbard U

    Result:

    it will generate results for a mu_list = np.arange(0.0, 3.1, 0.1)
    at fixed u.
    '''
    if os.path.isfile('result.h5'):
        return

    # set range of chemical potential mu.
    mu_list = np.arange(0.0, 3.1, 0.1)
    u_list = [u for mu in mu_list]
    generate_data(u_list, mu_list)


def plot_scan_u():
    with h5py.File('result.h5', 'r') as f:
        u_list = f['/u_list'][()]
        e_list = f['/e_list'][()]
        z_list = f['/z_list'][()]
        d_list = f['/d_list'][()]

    f, axarr = plt.subplots(3, sharex=True)
    axarr[0].plot(u_list, e_list)
    axarr[0].set_ylabel('energy')
    axarr[1].plot(u_list, d_list)
    axarr[1].set_ylabel('double occupancy')
    axarr[2].plot(u_list, z_list)
    axarr[2].set_ylabel('Z')
    axarr[2].set_xlabel('U')
    axarr[2].set_xlim(min(u_list), max(u_list))
    plt.show()
    f.savefig('result.png')


def plot_scan_mu():
    with h5py.File('result.h5', 'r') as f:
        mu_list = f['/mu_list'][()]
        e_list = f['/e_list'][()]
        z_list = f['/z_list'][()]
        d_list = f['/d_list'][()]
        n_list = f['/n_list'][()]

    f, axarr = plt.subplots(4, sharex=True)
    axarr[0].plot(mu_list, e_list)
    axarr[0].set_ylabel('energy')
    axarr[1].plot(mu_list, d_list)
    axarr[1].set_ylabel('double occupancy')
    axarr[2].plot(mu_list, z_list)
    axarr[2].set_ylabel('Z')
    axarr[3].plot(mu_list, n_list)
    axarr[3].set_ylabel('n')
    axarr[3].set_xlabel('$\mu$')
    axarr[3].set_xlim(min(mu_list), max(mu_list))
    plt.show()
    f.savefig('result.png')



if __name__=='__main__':
    if '-mu' in sys.argv:
        scan_mu()
        plot_scan_mu()
    else:
        scan_u()
        plot_scan_u()
