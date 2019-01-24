from __future__ import print_function
import numpy as np

try:
    from builtins import range, zip
except:
    pass


def fermi_dirac(e_fermi, delta, energy):
    """
    Return fermi-dirac distribution weight.
    """
    x = (energy - e_fermi)/delta
    if x < -200:
        f = 1.
    elif x > 200:
        f = 0.
    else:
        f = 1./(np.exp(x) + 1)
    return f


def num_electron_diff(e_fermi, delta, e_skn, w_k, nb_k, num_elec):
    ne = 0
    for e_kn in e_skn:
        for e_n, w, nb in zip(e_kn, w_k, nb_k):
            f = [fermi_dirac(e_fermi, delta, e) for e in e_n[:nb]]
            ne += np.sum(f)*w
    return ne - num_elec
