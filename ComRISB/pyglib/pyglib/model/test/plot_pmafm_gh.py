from scan_checkboard import get_scan_data
import matplotlib.pyplot as plt

def plot_scan_u_pmafm():
    '''plot a figure, which compares the PM and AFM results from Gutzwiller
    calculations.
    '''
    # Get data from precalcuated PM results.
    u_list_pm, e_list_pm, z_list_pm, d_list_pm, m_list_pm = \
            get_scan_data(fname='result_pm')

    # Get data from precalcuated AFM results.
    u_list_afm, e_list_afm, z_list_afm, d_list_afm, m_list_afm = \
            get_scan_data(fname='result_afm')

    # Get data from precalcuated Hartree-Fock AFM results.
    u_list_uhf, e_list_uhf, z_list_uhf, d_list_uhf, m_list_uhf = \
            get_scan_data(fname='result_afm_uhf')

    f, axarr = plt.subplots(2, 2, sharex=True)
    axarr[0, 0].plot(u_list_pm, e_list_pm, label='PM-G')
    axarr[0, 0].plot(u_list_afm, e_list_afm, label='AFM-G')
    axarr[0, 0].plot(u_list_uhf, e_list_uhf, label='AFM-HF')
    axarr[0, 0].legend()
    axarr[0, 0].set_ylabel('energy')
    axarr[1, 0].plot(u_list_pm, d_list_pm)
    axarr[1, 0].plot(u_list_afm, d_list_afm)
    axarr[1, 0].plot(u_list_uhf, d_list_uhf)
    axarr[1, 0].set_ylabel('double occupancy')
    axarr[0, 1].plot(u_list_pm, z_list_pm)
    axarr[0, 1].plot(u_list_afm, z_list_afm)
    axarr[0, 1].plot(u_list_uhf, z_list_uhf)
    axarr[0, 1].yaxis.tick_right()
    axarr[0, 1].yaxis.set_label_position("right")
    axarr[0, 1].set_ylabel('Z')
    axarr[1, 1].plot(u_list_pm, m_list_pm)
    axarr[1, 1].plot(u_list_afm, m_list_afm)
    axarr[1, 1].plot(u_list_uhf, m_list_uhf)
    axarr[1, 1].set_ylabel('local $<S_{z}>$')
    axarr[1, 1].set_ylim(-1,1)
    axarr[1, 1].yaxis.tick_right()
    axarr[1, 1].yaxis.set_label_position("right")
    axarr[1, 0].set_xlabel('U')
    axarr[1, 1].set_xlabel('U')
    axarr[1, 0].set_xlim(min(u_list_pm), max(u_list_pm))
    plt.tight_layout()
    plt.show()
    f.savefig('result_pmafm_gh.png')



if __name__=="__main__":
    plot_scan_u_pmafm()
