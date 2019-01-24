from builtins import zip
import h5py, numpy
import matplotlib.pyplot as plt
from pyglib.symm.unitary import comp_sph_harm_to_real_harm
from pyglib.basic.splot import colors
from matplotlib.gridspec import GridSpec

numpy.set_printoptions(precision=2, suppress=True)

u_csh2rh = comp_sph_harm_to_real_harm(orbital='d')
with h5py.File('GPARAM.h5', 'r') as f:
    u_csh2sab = f['/IMPURITY_1/DB_TO_SAB'][()].T
    u_csh2sab = u_csh2sab[:5, ::2]
u_rh2sab = u_csh2rh.T.conj().dot(u_csh2sab)

u_list = [0.0, 4.8, 4.8, 4.8, 4.8, 4.8, 4.8, 4.8, 4.8, 4.8, 4.8, 4.8]
j_list = [0.0, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

u_list = [4.8]
j_list = [0.8]

ef_list = []

for u, j in zip(u_list, j_list):
    with h5py.File('./U{}J{}/GBANDS_0.h5'.format(u, j), 'r') as f:
        nkpts = f['/IKP_END'][0]
        evals = []
        dwt_list = []

        with h5py.File('./U{}J{}/WH_RL_INP.h5'.format(u, j), 'r') as fz:
            R = fz['/IMPURITY_1/R'][::2, ::2].T

        for ik in range(nkpts):
            evals.append(f['/ISPIN_1/IKP_{}/ek'.format(ik+1)][()])
            evecs = f['/ISPIN_1/IKP_{}/ISYM_1/EVEC'.format(ik+1)][:, :5].T
            # Coherent part
            evecs = R.T.conj().dot(evecs)
            evecs = u_rh2sab.dot(evecs)
            dwt_list.append((evecs*evecs.conj()).real*20)
        evals = numpy.asarray(evals).T
        dwt_list = numpy.asarray(dwt_list)

    e_fermi = (numpy.max(evals[13]) + numpy.min(evals[14]))/2
    ef_list.append(e_fermi)


    print numpy.max(dwt_list)


    with h5py.File('../FeSb2Results.h5', 'r') as f:
        dos_e = f['/dos/u4.8j0.8/energies'][()]
        dos_t = f['/dos/u4.8j0.8/total'][()]
        dos_dxy = f['/dos/u4.8j0.8/dxy'][()]
        dos_dx2_y2 = f['/dos/u4.8j0.8/dx2-y2'][()]
        dos_dyz_dzx = f['/dos/u4.8j0.8/dyz+dzx'][()]
        dos_dz2 = f['/dos/u4.8j0.8/dz2'][()]


    evals -= e_fermi
    k_dist = numpy.loadtxt('wannier_band.dat').T[0]
    k_dist = k_dist[:nkpts]
    k_node = [0.00000, 0.99355, 1.53158, 2.52513, 3.06316, 3.54582, 4.08385,
            5.07740, 5.61543]
    k_label = ['$\Gamma$', 'Z', 'T', 'Y', '$\Gamma$', 'X', 'S', 'R', 'U']

    gs = GridSpec(1, 2, width_ratios=[2, 1])
    fig = plt.figure(figsize=(5,3))
    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1], sharey=ax1)

    for i, _eval in enumerate(evals):
        ax1.plot(k_dist, _eval, 'k-', linewidth=0.5)
        ax1.scatter(k_dist, _eval, s = dwt_list[:,1,i]/2
                +dwt_list[:,3,i]/2, c = 'b',
                edgecolors = 'b', alpha = 1)
        ax1.scatter(k_dist, _eval, s = dwt_list[:,0,i], c = 'r',
                edgecolors = 'r', alpha = 0.3)
    for kn in k_node:
        ax1.axvline(x=kn, linewidth=0.5, color='k', ls='--')
    ax1.axhline(y=0.0, linewidth=0.5, color='k', ls=':')
    ax1.set_ylim(-1.0, 1.0)
    ax1.set_xlim(k_dist[0], k_dist[-1])
    ax1.set_xticks(k_node)
    ax1.set_xticklabels(k_label)
    ax1.set_ylabel("Band energy (eV)")
    ax1.set_xlabel("k-path")

    ax2.get_yaxis().set_visible(False)
    ax2.fill_betweenx(dos_e, dos_t, 0., color=(0.7, 0.7, 0.7),
            facecolor=(0.7, 0.7, 0.7))
    ax2.plot(dos_dxy, dos_e, 'r-', label='$d_{xy}$')
    ax2.plot(dos_dyz_dzx, dos_e, 'b-', label='$d_{yz}+d_{zx}$')
    ax2.plot(dos_dx2_y2, dos_e, 'g-', label='$d_{x^{2}-y^{2}}$')
    ax2.plot(dos_dz2, dos_e, 'm-', label='$d_{z^{2}}$')
    ax2.legend(handlelength=1, handletextpad=0.2, fontsize='small')
    ax2.set_xlim(0, numpy.max(dos_t))
    ax2.set_xlabel('dos')
    fig.tight_layout()
    plt.subplots_adjust(wspace=0.05)
    plt.show()
    fig.savefig('bands_u{}j{}.pdf'.format(u,j))



print ef_list
