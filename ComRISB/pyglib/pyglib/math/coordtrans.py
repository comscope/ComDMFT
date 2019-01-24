import numpy as np


def cart2sph(v):
    '''converto v from cartesian coordinates to spherical coordinates.
    '''
    r = np.linalg.norm(v)
    abmax = np.max(np.abs(v[:2]))
    if abmax < 1.e-16:
        theta = 0.
        # take an arbitary value.
        phi = 0.
    else:
        v /= abmax
        ab = np.sqrt(v[0]**2+v[1]**2)
        phi = np.arcsin(v[1]/ab)
        if v[0] < 0.:
            phi = np.pi - phi
        theta = np.arccos(v[2]/np.linalg.norm(v))
    return r, theta, phi
