import numpy as np
import h5py
import pyglib.basic.units as units
import pyglib.basic.splot as splot

'''
Equation of state.
'''


def Murnaghan(parameters, vol):
    '''
    Given a vector of parameters and volumes, return a vector of energies.
    equation From PRB 28,5480 (1983)
    '''
    E0 = parameters[0]
    B0 = parameters[1]
    BP = parameters[2]
    V0 = parameters[3]
    return E0 + B0 * vol / BP * (((V0 / vol)**BP) /  \
            (BP - 1) + 1) - V0 * B0 / (BP - 1.0)


def Murnaghan_pv(parameters, vol):
    '''
    function P(V).
    '''
    B0 = parameters[1]
    BP = parameters[2]
    V0 = parameters[3]
    return B0 / BP * ((V0 / vol)**BP - 1.0)


def eos_fit_fun(pars, y, x):
    '''
    The objective function that will be minimized.
    '''
    return y - Murnaghan(pars, x)


def get_ev_fit(v, e):
    '''
    Fitting the Birch-Murnaghan EOS to data. v in \A^3, e in eV.
    Based on http://gilgamesh.cheme.cmu.edu/doc/software/jacapo/
    appendices/appendix-eos.html
    '''
    from pylab import polyfit
    from scipy.optimize import leastsq

    # fit a parabola to the data
    # y = ax^2 + bx + c
    a, b, c = polyfit(v, e, 2)
    '''The parabola does not fit the data very well, but we can use it to get
    some analytical guesses for other parameters.
    V0 = minimum energy volume, or where dE/dV=0
    E = aV^2 + bV + c
    dE/dV = 2aV + b = 0
    V0 = -b/2a
    E0 is the minimum energy, which is:
    E0 = aV0^2 + bV0 + c
    B is equal to V0*d^2E/dV^2, which is just 2a*V0
    and from experience we know Bprime_0 is usually a small number like 4
    '''
    # now here are our initial guesses.
    v0 = -b / (2 * a)
    e0 = a * v0**2 + b * v0 + c
    b0 = 2 * a * v0
    bP = 4
    # initial guesses in the same order used in the Murnaghan function
    x0 = [e0, b0, bP, v0]
    murnpars, ier = leastsq(eos_fit_fun, x0, args=(e, v))
    return murnpars


def h5get_mfit_ev(nmesh_fac=10, fsave='results.h5', path='/lapw'):
    '''Calculate and save Murnaghan fiting results in fsave.
    Interpolated e-v and p-v data on volume mesh with a factor a
    nmesh_fac of the original one are also stored.
    '''
    # Get e,v data.
    with h5py.File(fsave, 'r') as f:
        e_list = f[path+'/etot_list'][...]
        v_list = f['/vol_list'][...]

    # fitting
    murnpars = get_ev_fit(v_list, e_list)

    vh = np.linspace(v_list[0], v_list[-1], nmesh_fac * len(v_list) - 1)
    eh = Murnaghan(murnpars, vh)
    ph = Murnaghan_pv(murnpars, vh)*units.eVA_GPa
    with h5py.File(fsave, 'a') as f:
        if path+'/eosfit' in f:
            del f[path+'/eosfit']
        f[path+'/eosfit/e0'] = murnpars[0]
        f[path+'/eosfit/b0'] = murnpars[1]
        f[path+'/eosfit/bp'] = murnpars[2]
        f[path+'/eosfit/v0'] = murnpars[3]
        f[path+'/eosfit/v_list'] = vh
        f[path+'/eosfit/e_list'] = eh
        f[path+'/eosfit/p_list'] = ph
    splot.xy2_plot([v_list, vh], [e_list, eh], ['o', '-'], ['raw', 'fitting'],
            xlabel='V ($\AA^3$/primitive cell)',
            ylabel='E (eV/primitive cell)', fsave=path+'_evfit.pdf')
    splot.xy_plot(vh, ph, xlabel='V ($\AA^3$/primitive cell)',
            ylabel='P (GPa)', fsave=path+'_pvfit.pdf')


def eos_spline(v, e, tol):
    '''
    Get volume, energy, pressure, and bulk modulus using spline, given
    v in \A^3 and e in eV.
    '''
    from scipy.interpolate import UnivariateSpline
    s = UnivariateSpline(v, e, k=3, s=tol)
    vh = np.linspace(v[0], v[-1], 10 * len(v) - 1)
    eh = [s.derivatives(i)[0] for i in vh]
    ph = [-s.derivatives(i)[1] * units.eVA_GPa for i in vh]
    bh = [s.derivatives(i)[2] * vh[i] * units.eVA_GPa for i in vh]
    return vh, eh, ph, bh
