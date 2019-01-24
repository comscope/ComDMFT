import numpy as np
from sympy.physics.quantum.cg import CG


def jj_to_cubic_relativistic_harmonics(orbital='f'):
    if 'f' == orbital:
        jj_to_cubic = np.zeros((14,14))
        jj_to_cubic[0,8] = -np.sqrt(1./6.) # |5/2, -5/2>
        jj_to_cubic[4,8] =  np.sqrt(5./6.) # |5/2, +3/2> G7, 5/2, +
        jj_to_cubic[5,10] = -np.sqrt(1./6.) # |5/2, +5/2>
        jj_to_cubic[1,10] =  np.sqrt(5./6.) # |5/2, -3/2> G7, 5/2, -

        jj_to_cubic[4,0] =  np.sqrt(1./6.) # |5/2, +3/2>
        jj_to_cubic[0,0] =  np.sqrt(5./6.) # |5/2, -5/2> G81, 5/2, +
        jj_to_cubic[1,2] =  np.sqrt(1./6.) # |5/2, -3/2>
        jj_to_cubic[5,2] =  np.sqrt(5./6.) # |5/2, +5/2> G81, 5/2, -

        jj_to_cubic[3,4] =  1. # |5/2, +1/2> G82, 5/2, +
        jj_to_cubic[2,6] =  1. # |5/2, -1/2> G82, 5/2, -

        jj_to_cubic[13,12] =  np.sqrt(5./12.) # |7/2, +7/2>
        jj_to_cubic[ 9,12] =  np.sqrt(7./12.) # |7/2, -1/2> G6, 7/2, +
        jj_to_cubic[ 6,13] =  np.sqrt(5./12.) # |7/2, -7/2>
        jj_to_cubic[10,13] =  np.sqrt(7./12.) # |7/2, +1/2> G6, 7/2, -

        jj_to_cubic[12,11] = -np.sqrt(3./4.) # |7/2, +5/2>
        jj_to_cubic[ 8,11] =  np.sqrt(1./4.) # |7/2, -3/2> G7, 7/2, +
        jj_to_cubic[ 7,9] =  np.sqrt(3./4.) # |7/2, -5/2>
        jj_to_cubic[11,9] = -np.sqrt(1./4.) # |7/2, +3/2> G7, 7/2, -

        jj_to_cubic[13,7] =  np.sqrt(7./12.) # |7/2, +7/2>
        jj_to_cubic[ 9,7] = -np.sqrt(5./12.) # |7/2, -1/2> G81, 7/2, +
        jj_to_cubic[ 6,5] = -np.sqrt(7./12.) # |7/2, -7/2>
        jj_to_cubic[10,5] =  np.sqrt(5./12.) # |7/2, +1/2> G81, 7/2, -

        jj_to_cubic[12,3] = -np.sqrt(1./4.)  # |7/2, +5/2>
        jj_to_cubic[ 8,3] = -np.sqrt(3./4.)  # |7/2, -3/2> G82, 7/2, +
        jj_to_cubic[ 7,1] =  np.sqrt(1./4.)  # |7/2, -5/2>
        jj_to_cubic[11,1] =  np.sqrt(3./4.)  # |7/2, +3/2> G82, 7/2, -
    else:
        raise ValueError('UndefinedFunction')
    return jj_to_cubic


def comp_sph_harm_to_real_harm(dim_m):
    csh2rh = np.zeros((dim_m, dim_m), dtype=complex)
    # qunatum number l
    l = dim_m//2
    # set the orbital index
    iy = range(dim_m)
    for i in range(l):
        iy += [iy.pop(0)]
    for m in range(-l,l+1):
        if m < 0:
            csh2rh[iy[m], iy[m]] = 1.j/np.sqrt(2.)
            csh2rh[iy[-m], iy[m]] = -1.j/np.sqrt(2.)*(-1)**m
        elif m == 0:
            csh2rh[iy[m], iy[m]] = 1.
        else:
            csh2rh[iy[-m], iy[m]] = 1./np.sqrt(2)
            csh2rh[iy[m], iy[m]] = 1./np.sqrt(2)*(-1)**m
    return csh2rh


def get_u_csh2rh_all(ncorbs_list):
    u_csh2rh_list = [comp_sph_harm_to_real_harm(ncorbs) for \
            ncorbs in ncorbs_list]
    return u_csh2rh_list


def comp_sph_harm_to_relativistic_harm(dim_ms):
    '''transformation matrix from spin-complex spherical harmonics
    (orbital fast) to relativistic harmonics
    '''
    csh2relh = np.zeros((dim_ms, dim_ms), dtype=complex)
    dim_m = dim_ms//2
    l = dim_m//2
    # set the orbital index
    iy = range(dim_m)
    for i in range(l):
        iy += [iy.pop(0)]
    # add slow spin index
    # for spin: up, then dn. wien2k convention.
    iys = {0.5:iy, -0.5:[iy[i]+dim_m for i in range(dim_m)]}
    # j=l-1/2 block
    # relativistic_harmonics index
    i_jm = -1
    for i in [-0.5,0.5]:
        _j = l+i
        for mj in np.arange(-_j, _j+1):
            i_jm += 1
            for s in [-0.5,0.5]:
                csh2relh[iys[s][int(round(mj-s))], i_jm] = \
                        CG(l,mj-s,0.5,s,_j,mj).doit()
    return csh2relh


def get_u_csh2relh_all(ncorbs_list):
    u_csh2relh_list = [comp_sph_harm_to_relativistic_harm(ncorbs) for \
            ncorbs in ncorbs_list]
    return u_csh2relh_list


def get_u_csh2wan_all(ncorbs_list):
    ncorbs = ncorbs_list[0]
    if ncorbs%2 == 0: #with spin-orbit interaction
        u_csh2wan_list = get_u_csh2relh_all(ncorbs_list)
    else:
        u_csh2wan_list = get_u_csh2rh_all(ncorbs_list)
    return u_csh2wan_list
