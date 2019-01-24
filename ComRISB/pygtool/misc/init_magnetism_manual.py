import numpy
from pyglib.gutz.init_magnet_conf import init_magnet_conf_with_init_dm

# create a list of initial density matrices
dm_list = []

# the first impurity is Ce with spin-orbit dim = 7*2 = 14.
dm = numpy.zeros((14,14))

# put 1 f-electron in |J=5/2, Jz=-5/2>.
dm[0, 0] = 1.
dm_list.append(dm)

# generate the initial input file `GVEXT.h5`
init_magnet_conf_with_init_dm(dm_list)

