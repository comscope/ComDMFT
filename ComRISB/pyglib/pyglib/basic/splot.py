from __future__ import absolute_import, division, print_function
try:
    from builtins import (bytes, str, open, super, range,
            zip, round, input, int, pow, object)
except:
    pass

import matplotlib.pyplot as plt


'''Simple pdf plot for quick visual check.
'''

colors = ['black', 'red', 'green', 'blue', 'orange','violet',
        'darkred', 'darkgreen', 'navy', 'brown', 'chocolate', 'darkorange',
        'gold', 'olive', 'maroon']

def xy_plot(x_list, y_list, xlabel='x', ylabel='y', fsave='test.pdf'):
    plt.figure()
    plt.plot(x_list, y_list)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.savefig(fsave)


def xy2_plot(x_list, y_list, pattern_list, label_list,
        xlabel='x', ylabel='y', fsave='test.pdf'):
    plt.figure()
    for x, y, pattern, label in zip(x_list, y_list, pattern_list, label_list):
        plt.plot(x, y, pattern, label=label)
    plt.legend()
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.savefig(fsave)
