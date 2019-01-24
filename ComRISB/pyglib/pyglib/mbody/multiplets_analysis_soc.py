################################################
#Multiplets analysis for LDA+G+SOC calculations.
################################################

import numpy as np
import h5py
from pyglib.mbody.multiplets_analysis_lib import get_local_histogram


def calc_save_atomic_states(imp=1, num_ev=100):
    r'''
    Calculate and save the eigen-values of the
    local many-body density matrix
    and the labels in file 'multiplets.h5'.

    Parameters:

    * imp: int
        one-based index of impurity
    * num_ev: int
        number of biggest eigenvalues to be calculated

    Note: Results stored in 'multiplets.h5' file, with entries:

    * sum_weights: float
        sum of the list of weight calculated.
    * weights: float array
        array of the weight
    * n_labels: int array
        array of electron occupancy
    * s_labels: float array
        array of s spin index
    * l_labels: float array
        array of l angular momentum index
    * j_labels: float array
        array of j total angular momentum index
    * eval_s: float
        expectation value of s operator
    * eval_l: float
        expectation value of l operator
    * eval_j: float
        expectation value of j operator
    * ent_entropy: float
        entanglement entropy
    * deg_labels: int array
        array of degeneracy.
    '''
    # Get eigen-values and labels of the local atomic states.
    vals, n_labels, s_labels, l_labels, j_labels, \
            chi_labels, multiplet_degeneracies, \
            eval_s, eval_l, eval_j, ent_entropy \
            = get_local_histogram(imp, num_ev=num_ev)

    # Store results to metadata.
    with h5py.File("multiplets.h5".format(imp), 'a') as f:
        base = "/impurity_" + str(imp)
        if base in f:
            del f[base]
        f[base + "/sum_weights"] = np.sum(vals)
        f[base + "/weights"] = vals
        f[base + "/n_labels"] = n_labels
        f[base + "/s_labels"] = s_labels
        f[base + "/l_labels"] = l_labels
        f[base + "/j_labels"] = j_labels
        f[base + "/eval_s"] = eval_s
        f[base + "/eval_l"] = eval_l
        f[base + "/eval_j"] = eval_j
        f[base + "/ent_entropy"] = ent_entropy
        f[base + "/deg_labels"] = multiplet_degeneracies


def plot_atomic_states(imp=1, num_label=5):
    r'''
    Plot the atomic state probabilities using the formula
    :math:`\rho=e^{-F}`.

    Parameters

    * imp: int
        one-based impurity index
    * num_label: int
        number of significant configurations to be labeled in the figure.

    Note: A figure will be saved in 'histogram.pdf' file.
    '''
    # Read in metadata
    with h5py.File("multiplets.h5", 'r') as f:
        base = "/impurity_" + str(imp)
        vals = f[base + "/weights"][...]
        multiplet_degeneracies = f[base + "/deg_labels"][...]
        n_labels = f[base + "/n_labels"][...]
        j_labels = f[base + "/j_labels"][...]

    # Generate data corresponding to \rho = exp(-F)
    f_wt = -np.log(vals / multiplet_degeneracies)
    f_wt = f_wt - f_wt[0]

    # Labels of the atomic states
    labels = ['N=%1d,J=%3.1f' % (int(n + 0.1), j) for n, j in
              zip(n_labels[:num_label], j_labels[:num_label])]

    # Scatter plot with annotation.
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    ax.scatter(f_wt, vals)
    for label, x, y in zip(labels, f_wt[:num_label], vals[:num_label]):
        ax.annotate(
            label,
            xy=(x, y), xytext=(0, 20),
            textcoords='offset points', ha='center', va='bottom',
            bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.5),
            arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0'))
    ax.set_xlim(-0.1, 12)
    ax.set_ylim(0, 1)
    ax.set_xlabel("$f_{n}$")
    ax.set_ylabel("$e^{-f_{n}}d_{n}/\sum_{n}{e^{-f_{n}}d_{n}}$")
    plt.title("Eigen-values of the local many-body density matrix")
    fig.tight_layout()
    plt.show()
    fig.savefig('histogram.pdf')


if __name__ == "__main__":
    calc_save_atomic_states(imp=1, num_ev=100)
    plot_atomic_states(imp=1, num_label=5)
