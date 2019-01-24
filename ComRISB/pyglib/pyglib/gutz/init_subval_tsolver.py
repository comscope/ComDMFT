#!/usr/bin/env python
from __future__ import print_function
import numpy,h5py
from pyglib.gutz.usrqa import get_usr_input


def init_subval_tsolver():
    '''setup the parameters for sub-valence truncation solver.
    '''
    with h5py.File('GPARAM.h5', 'r') as f:
        imap_imp = f['/IMAP_IMP'][()]
        gpms_list = []
        ngpm_list = []
        nval_list = []
        ngp_list = []
        for i, imap in enumerate(imap_imp):
            # one-base fortran index
            imp = i+1
            if imap < imp:
                gpms_list.append(gpms_list[imap])
                ngpm_list.append(ngpm_list[imap])
                nval_list.append(nval_list[imap])
                ngp_list.append(ngp_list[imap])
            elif imap == imp:
                sx = f['/IMPURITY_{}/SX'.format(imp)][()].T
                sy = f['/IMPURITY_{}/SY'.format(imp)][()].T
                sz = f['/IMPURITY_{}/SZ'.format(imp)][()].T
                lx = f['/IMPURITY_{}/LX'.format(imp)][()].T
                ly = f['/IMPURITY_{}/LY'.format(imp)][()].T
                lz = f['/IMPURITY_{}/LZ'.format(imp)][()].T
                jx, jy, jz = sx + lx, sy + ly, sz + lz
                j2 = jx.dot(jx) + jy.dot(jy) + jz.dot(jz)
                # assuming 1-l channel
                l = (sx.shape[0]/2-1)/2
                ngpm = numpy.zeros([2], dtype=numpy.int)
                gpms = numpy.zeros([2, 2*l+2], dtype=numpy.int)
                j21, j22 = (l-0.5)*(l+0.5), (l+0.5)*(l+1.5)
                for ii, val in enumerate(j2.diagonal()):
                    if abs(val - j21) < 1.e-6:
                        # for frotran, one-base
                        gpms[0, ngpm[0]] = ii + 1
                        ngpm[0] += 1
                    elif abs(val - j22) < 1.e-6:
                        # for frotran, one-base
                        gpms[1, ngpm[1]] = ii + 1
                        ngpm[1] += 1
                    else:
                        raise ValueError( \
                                ' inconsistent setup with j2-diagonal = {}'\
                                .format(val))
                print(' IMPURITY {}\n num orbs groups = {}'.format(imp, 2))
                nvals = []
                for ig in range(2):
                    print(' group {} member indices: {}'.format(ig, \
                            gpms[ig, :ngpm[ig]]))
                    nval_bot = get_usr_input( \
                            ' choose min. occupancy of physical subspace', \
                            map(str, range(ngpm[ig]+1)))
                    nval_bot = int(nval_bot)
                    nval_top = get_usr_input( \
                            ' choose min. occupancy of physical subspace', \
                            map(str, range(nval_bot, ngpm[ig]+1)))
                    nval_top = int(nval_top)
                    nvals.append([nval_bot, nval_top])
                    nval_bot = get_usr_input( \
                            ' choose min. occupancy of variational subspace', \
                            map(str, range(ngpm[ig]+1)))
                    nval_bot = int(nval_bot)
                    nval_top = get_usr_input( \
                            ' choose min. occupancy of variational subspace', \
                            map(str, range(nval_bot, ngpm[ig]+1)))
                    nval_top = int(nval_top)
                    nvals.append([nval_bot, nval_top])
                nval_list.append(nvals)
                ngpm_list.append(ngpm)
                gpms_list.append(gpms)
                ngp_list.append(2)
            else:
                raise ValueError(' imap = {} > imp = {}'.format(imap, imp))

    with h5py.File('GESOLVER.h5', 'w') as f:
        for i, ngp in enumerate(ngp_list):
            imp = i + 1
            f['/IMPURITY_{}/n_orbs_group'.format(imp)] = [ngp]
            f['/IMPURITY_{}/n_gpmembers'.format(imp)] = ngpm_list[i]
            f['/IMPURITY_{}/gpmember_list'.format(imp)] = gpms_list[i]
            f['/IMPURITY_{}/gpnocc_range'.format(imp)] = nval_list[i]


if __name__ == '__main__':
    init_subval_tsolver()
