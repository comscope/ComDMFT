from __future__ import print_function
'''Script works for both metallic phase and Mott phase with
spin-orbit interaction ONLY.'''

import os, sys
import subprocess as sp


msg = r'''local multiplet analysis for both metallic phase and mott phase.
    possible labels includ N, S, L and J, dependent of the specific symmetry.
    inline command arguments:
    -i int -- specify the one-based impurity index
    -f str -- specify the hdf5 file name storing multiplet analysis results.
'''

if '-h' in sys.argv:
    print(msg)
    sys.exit()

# impurity index
if '-i' in sys.argv:
    imp = int(sys.argv[sys.argv.index('-i')+1])
else:
    imp = 1

# hdf5 file name of the multiplet analysis results.
if '-f' in sys.argv:
    fname = sys.argv[sys.argv.index('-f')+1]
else:
    fname = 'multiplets.h5'

if not os.path.isfile('EMBED_HAMIL_ANALYSIS_{}.h5'.format(imp)):
    cmd = [os.environ['WIEN_GUTZ_ROOT2']+'/exe_spci_j2_mott_analysis', \
            str(imp)]
    sp.call(cmd)

from pyglib.mbody.multiplets_analysis_s_square import calc_save_atomic_states, \
        plot_atomic_states

# Calculate local histograms for impurity 1
# Check the 100 dominant eigen-states.
calc_save_atomic_states(imp=imp, num_ev=400, fname=fname)

# Plot the histogram
plot_atomic_states(imp=imp, num_label=7)
