from __future__ import print_function
import numpy
from builtins import zip
from scipy.special import erfc
from scipy.optimize import bisect


def get_fermi_level(bnd_es, wklist, num_e, delta=0.0258, \
        ismear=0, iso=1):
    # delta = 0.0258519909eV = 300K
    emin = numpy.min(bnd_es)
    emax = numpy.max(bnd_es)
    efermi = bisect(err_fun, emin, emax, \
            args=(bnd_es, wklist, num_e, delta, ismear, iso))
    return efermi


def err_fun(mu, bnd_es, wklist, num_e, delta, ismear, iso):
    ferwes = get_fermi_weight(mu, bnd_es, wklist, delta, ismear=ismear,
            iso=iso)
    _nume = numpy.sum(ferwes)
    return _nume - num_e


def get_fermi_weight(mu, bnd_es, wklist, delta=0.0258, \
        ismear=0, iso=1):
    ferwes = []
    for bnd_e in bnd_es:
        ferwes.append([])
        for bnd_ek, wk in zip(bnd_e, wklist):
            ferwes[-1].append([])
            for e in bnd_ek:
                x = (e-mu)/delta
                if ismear == -1:
                    fw = fermi_dirac(x)
                elif ismear == 0:
                    fw = gaussian(x)
                else:
                    raise ValueError("Not defined ismear = {}!".format(ismear))
                fw *= wk
                ferwes[-1][-1].append(fw)
    ferwes = numpy.asarray(ferwes)
    if len(ferwes) == 1 and iso == 1:
        ferwes *= 2
    return ferwes


def fermi_dirac(x):
    """
    Return fermi-dirac distribution weight.
    """
    if x < -200:
        f = 1.
    elif x > 200:
        f = 0.
    else:
        f = 1./(numpy.exp(x) + 1)
    return f


def gaussian(x):
    """
    Return gaussian distribution weight.
    """
    if x < -7:
        f = 2.
    elif x > 7:
        f = 0.
    else:
        f = erfc(x)
    return f/2.0
