from __future__ import print_function
import matplotlib.pyplot as plt
try:
    from builtins import zip
except:
    pass

import numpy, h5py
from pyglib.basic.splot import colors


def s_to_index(s):
    '''get integer index for s, which can either be integer or halves.
    '''
    if 0.49 < s % 1 < 0.51: # half s, 1/2, 3/2, ...
        return int(s - 0.49)
    elif 0 <= s % 1 < 0.01 or 0.99 < s % 1: # integer s, 0, 1, ...
        return int(s + 0.1)
    else:
        raise ValueError(' s = {}'.format(s))


def get_group_weights(weights_list, n_labels_list, s_labels_list):
    '''group weights according to valence n and s quantum number.
    '''
    nmin = numpy.min(n_labels_list[0])
    nmax = numpy.max(n_labels_list[0])
    smax = numpy.max(s_labels_list[0])
    ns = s_to_index(smax)

    g_weights_list = []
    for weights, n_labels, s_labels in zip(weights_list, n_labels_list,
            s_labels_list):
        g_weights_list.append(numpy.zeros([nmax-nmin+1, ns+1]))
        for wt, n, s in zip(weights, n_labels, s_labels):
            if wt < 1.e-7:
                continue
            try:
                i = s_to_index(s)
            except ValueError:
                print(' wt = {}: to be ignored!'.format(wt))
                continue
            g_weights_list[-1][n-nmin, i] += wt
    return g_weights_list


def hist_ns(weights_list, n_labels_list, s_labels_list,
        amlabel='S', tol=1.e-2):
    '''generate histogram plot for the multiple local reduced density matrix
    with label n and s(j).
    '''
    x_ticks = []
    label_ticks = []
    icolor = 0
    color_map = [{}, {}]
    nmin = numpy.min(n_labels_list[0])
    g_weights_list = get_group_weights(weights_list, n_labels_list,
            s_labels_list)
    idx = -1
    fig, ax = plt.subplots(figsize=(6,3))
    g_weights_list = numpy.asarray(g_weights_list)
    rwidth = 1./(1+g_weights_list.shape[0])

    cmap = plt.get_cmap('plasma')

    nalpha_start = -1
    for n in range(g_weights_list.shape[1]):
        wmax = numpy.max(g_weights_list[:, n, :])
        if wmax >= tol:
            if nalpha_start == -1:
                nalpha_start = n
            else:
                nalpha_end = n

    for n in range(g_weights_list.shape[1]):
        idx_start = idx
        for s in range(g_weights_list.shape[2]):
            wmax = numpy.max(g_weights_list[:, n, s])
            jn = (n+nmin) % 2
            if wmax >= tol:
                idx += 1
                if s in color_map[jn]:
                    color = color_map[jn][s]
                    label = None
                else:
                    icolor += 1
                    color = color_map[jn][s] = colors[icolor]
                    if jn == 0:
                        label = '{}={}'.format(amlabel, s)
                    else:
                        label = '{}={}/2'.format(amlabel, 2*s+1)
                for i, w in enumerate(g_weights_list[:, n, s]):
                    ax.hist([idx-0.5+rwidth*(i+1)], 1, weights=[w],
                            rwidth=rwidth*0.9, align='mid',
                            color=color, label=label)
                    label = None

        if idx > idx_start:
            ax.axvspan(idx_start+0.5, idx+0.5,
                    facecolor=cmap(float(n-nalpha_start) \
                    /(nalpha_end-nalpha_start)),
                    zorder=-1, alpha=0.2, ls='--', ec='blue')
            x_ticks.append((idx_start+idx+1)/2.)
            label_ticks.append('N='+str(n+nmin))
    ax.set_xlim(-0.5, idx+0.5)
    ax.set_xticks(x_ticks)
    ax.set_xticklabels(label_ticks)
    ax.legend()
    ax.set_ylabel('Probability')
    fig.tight_layout()
    plt.show()
    fig.savefig('hist_s.pdf')



if __name__=='__main__':
    weights_list = []
    n_labels_list = []
    s_labels_list = []
    for i in range(2):
        with h5py.File('multiplets_{}.h5'.format(i), 'r') as f:
            weights, n_labels, s_labels = f['/impurity_1/weights'][()], \
                    f['/impurity_1/n_labels'][()], \
                    f['/impurity_1/s_labels'][()]
            weights_list.append(weights)
            n_labels_list.append(n_labels)
            s_labels_list.append(s_labels)
    hist_ns(weights_list, n_labels_list, s_labels_list)
