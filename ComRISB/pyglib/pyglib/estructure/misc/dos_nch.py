from math import pi, sqrt
import numpy as np
try:
    from builtins import zip
except:
    pass

from itertools import product,count
from pyglib.estructure import fermi
from pyglib.basic import units


class DOS:

    def __init__(self, w_k, e_skn,  width=0.1, window=None, npts=201):
        """Electronic Density Of States object.

        width: float
          Width of guassian smearing.
        window: tuple of two float
          Use ``window=(emin, emax)``.  If not specified, a window
          big enough to hold all the eigenvalues will be used.
        npts: int
          Number of points.

        """

        self.npts = npts
        self.width = width
        self.w_k = w_k
        self.e_skn = e_skn

        if window is None:
            emin = self.e_skn.min() - 5 * self.width
            emax = self.e_skn.max() + 5 * self.width
        else:
            emin, emax = window

        self.energies = np.linspace(emin, emax, npts)


    def delta(self, energy):
        """
        Return a delta-function centered at energy.
        """
        x = -((self.energies - energy) / self.width)**2
        return np.exp(x) / (sqrt(pi) * self.width)


    def get_energies(self):
        """"
        return energy mesh.
        """
        return self.energies


    def get_dos_component(self, psi_skn_w):
        """
        Get array of DOS values.
        """
        dos_list = []

        # spin
        for e_kn, psi_kn_w in zip(self.e_skn,psi_skn_w):
            dos = np.zeros(self.npts, dtype=np.complex)
            # k-point
            for w, e_n, psi_n_w in zip(self.w_k, e_kn, psi_kn_w):
                # band index
                for e, psi_w in zip(e_n, psi_n_w):
                    dos += w * self.delta(e) * psi_w
            dos_list.append(dos)
        return np.array(dos_list)


    def get_dos_t(self):
        """
        Get array of DOS values.
        """
        dos_list = []

        # spin
        for e_kn in self.e_skn:
            dos = np.zeros(self.npts)
            # k-point
            for w, e_n in zip(self.w_k, e_kn):
                # band index
                for e in e_n:
                    dos += w * self.delta(e)
            dos_list.append(dos)
        return np.array(dos_list)


    def get_dos_noncoherent(self, psi_sksn, ecdfdffdelta, f_skn, bnd_ne,
            i):
        """
        Non-coherent part of dos.
        """

        # symmetry operations
        s_range = range(psi_sksn.shape[2])
        norb = ecdfdffdelta.shape[0]

        dos_list = []
        # spin
        for e_kn, f_kn, psi_ksn in zip(self.e_skn, f_skn, psi_sksn):
            dos = np.zeros(self.npts, dtype=np.complex)

            # k-point
            for e_n1, f_n1, psi_sn1, nb1 in zip(e_kn, f_kn, psi_ksn, bnd_ne):
                f_k1 = f_n1[0]
                for n1, e1, f1 in zip(count(),
                        e_n1[nb1[1]:nb1[2]+1],
                        f_n1[nb1[1]:nb1[2]+1]):
                    if(abs(e1) > 10): continue
                    for n2, e2, f2 in zip(count(),
                            e_n1[nb1[1]:nb1[2]+1],
                            f_n1[nb1[1]:nb1[2]+1]):
                        if(abs(e2) > 10): continue
                        for n3, e3, f3 in zip(count(),
                                e_n1[nb1[1]:nb1[2]+1],
                                f_n1[nb1[1]:nb1[2]+1]):
                            if(abs(e3) > 10): continue
                            res1 = (f_k1-f1)*(f_k1-f2)*f3 + f1*f2*(f_k1-f3)
                            if(abs(res1)<1.e-16):
                                continue
                            if(abs(f1/f_k1 -0.5) > 0.0001 and
                                    abs(f2/f_k1 -0.5) > 0.0001 and
                                    abs(f3/f_k1 -0.5) > 0.0001):
                                continue

                            res1 *= self.delta(e1+e2-e3)
                            for s1 in s_range:
                                resi = np.einsum('abc,b,c,a',
                                        ecdfdffdelta[i],
                                        psi_sn1[s1,n1,:norb]\
                                        .conj(),
                                        psi_sn1[s1,n2,:norb] \
                                        .conj(),
                                        psi_sn1[s1,n3,:norb])
                                dos += res1*np.conj(resi)*resi
            dos_list.append(dos)
        return np.array(dos_list)


    def get_dos_noncoherent_all(self, psi_sksn, ecdfdffdelta, f_skn, bnd_ne,
            i, j):
        """
        Non-coherent part of dos.
        """

        # symmetry operations
        s_range = range(psi_sksn.shape[2])
        norb = ecdfdffdelta.shape[0]

        dos_list = []
        # spin
        for e_kn, f_kn, psi_ksn in zip(self.e_skn, f_skn, psi_sksn):
            dos = np.zeros(self.npts, dtype=np.complex)

            # k-point
            for e_n1,f_n1, psi_sn1, nb1 in zip(e_kn, f_kn, psi_ksn, bnd_ne):
                # k-point
                for e_n2, f_n2, psi_sn2, nb2 in zip(e_kn, f_kn, psi_ksn,
                        bnd_ne):
                    # k-point
                    for e_n3, f_n3, psi_sn3, nb3 in zip(e_kn, f_kn, psi_ksn,
                            bnd_ne):
                        for n1, e1, f1 in zip(count(),
                                e_n1[nb1[1]:nb1[2]+1],
                                f_n1[nb1[1]:nb1[2]+1]):
                            for n2, e2, f2 in zip(count(),
                                    e_n2[nb2[1]:nb2[2]+1],
                                    f_n2[nb2[1]:nb2[2]+1]):
                                    for n3, e3, f3 in zip(count(),
                                            e_n3[nb3[1]:nb3[2]+1],
                                            f_n3[nb3[1]:nb3[2]+1]):
                                        res1 = (1-f1)*(1-f2)*f3 + f1*f2*(1-f3)
                                        res1 *= self.delta(e1+e2-e3)
                                        for s1, s2, s3 in product(s_range,
                                                s_range,s_range):
                                            resi = np.einsum('abc,b,c,a',
                                                    ecdfdffdelta[i],
                                                    psi_sn1[s1,n1,:norb]\
                                                    .conj(),
                                                    psi_sn2[s2,n2,:norb] \
                                                    .conj(),
                                                    psi_sn3[s3,n2,:norb])
                                            resj = np.einsum('abc,b,c,a',
                                                    ecdfdffdelta[j],
                                                    psi_sn1[s1,n1,:norb] \
                                                    .conj(),
                                                    psi_sn2[s2,n2,:norb] \
                                                    .conj(),
                                                    psi_sn3[s3,n3,:norb])
                                            dos += res1*np.conj(resi)*resj
            dos_list.append(dos)
        return np.array(dos_list)


def get_all_psi_skn_w(e_skn, psi_sksn, bnd_ne):
    '''
    Reduce psi_sksn to psi_skn_w.
    '''
    psi_skn_w = np.zeros(list(e_skn.shape) + [psi_sksn.shape[-1]], dtype=float)
    for isp in range(e_skn.shape[0]):
        for k in range(e_skn.shape[1]):
            for isy in range(psi_sksn.shape[2]):
                for n in range(bnd_ne[k, 1], bnd_ne[k, 2]):
                    for a in range(psi_sksn.shape[-1]):
                        psi_skn_w[isp, k, n, a] += (
                                psi_sksn[isp, k, isy, n - bnd_ne[k, 1], a] *
                                np.conj(
                                psi_sksn[isp, k, isy, n - bnd_ne[k, 1], a])).\
                                real
    return psi_skn_w / psi_sksn.shape[2]


def get_all_psi_skn_w_ab(e_skn, psi_sksn, bnd_ne):
    '''
    Reduce psi_sksn to psi_skn_w_ab.
    '''
    n_orb = psi_sksn.shape[-1]
    psi_skn_w_ab = np.zeros(list(e_skn.shape) + [n_orb, n_orb], dtype=complex)
    for k in range(e_skn.shape[1]):
        for n in range(bnd_ne[k, 1], bnd_ne[k, 2]):
            psi_skn_w_ab[:, k, n, :, :] = np.einsum('ijk, ijl->ikl',
                    psi_sksn[:, k, :, n - bnd_ne[k, 1], :],
                    np.conj(
                    psi_sksn[:, k, :, n - bnd_ne[k, 1], :]))
    return psi_skn_w_ab / psi_sksn.shape[2]



if __name__ == "__main__":
    '''
    Test run.
    '''
    import h5py

    with h5py.File("GPARAMBANDS.h5", 'r') as f:
        w_k = f["/kptwt"][...]
        bnd_ne = f["/NE_LIST"][...]
        nsym = f["/symnop"][...]
        nkp = len(w_k)

    with h5py.File("GLOG.h5", 'r') as f:
        e_fermi = f["/e_fermi"][0]
        nspin = f["/nspin"][0]

    with h5py.File("GBANDS_0.h5", 'r') as f:
        e_skn = []
        for isp in range(nspin):
            e_kn = []
            for k in range(nkp):
                e_kn.append(f['/ISPIN_{}/IKP_{}/ek'.format(isp+1,k+1)][...])
            e_skn.append(e_kn)

    num_elec = np.sum(f_skn) - 2
    print(fermi.num_electron_diff(0.0, 0.001, e_skn, w_k, bnd_ne[:,0], num_elec))

    # test fermi level
    from scipy.optimize import newton
    delta_list = [10.0, 8.0, 6.0, 4.0, 3.0, 2.0, 1.0, 0.9, 0.8, 0.6, 0.5, 0.3, 0.2, 0.1, 0.01]
    f_list = []
    x0 = 0
    for delta in delta_list:
        x0 = newton(fermi.num_electron_diff, x0, \
                args = (delta, e_skn, w_k, bnd_ne[:,0], num_elec))
        print x0

    psi_skn_w_ab = get_all_psi_skn_w_ab(e_skn, psi_sksn, bnd_ne)
    window = (-5., 15.)
    dos = DOS(w_k, e_skn,  width=0.05, window=window, npts=1001)
    energies = dos.get_energies()
    dos_t = dos.get_dos_t()

    # dos_f
    psi_sksn_f = np.einsum('...ii', psi_skn_w_ab[...,:,:])
    dos_f = dos.get_dos_component(psi_sksn_f)

    import matplotlib.pyplot as plt

    plt.figure(0)
    plt.plot(energies, dos_t[0].real, '-', label='t')
    plt.plot(energies, dos_f[0].real, 'o', label='f')
    plt.legend()
    plt.xlabel("E (eV)")
    plt.ylabel("DOS (states/f.u.)")
#    plt.savefig('dos.pdf')
    plt.show()
