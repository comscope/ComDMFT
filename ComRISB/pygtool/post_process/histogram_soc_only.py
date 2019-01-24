'''Script works for both metallic phase and Mott phase with
spin-orbit interaction ONLY.'''

import os
import subprocess as sp

# impurity index
imp = 1

if not os.path.isfile('EMBED_HAMIL_ANALYSIS_{}.h5'.format(imp)):
    cmd = [os.environ['WIEN_GUTZ_ROOT2']+'/exe_spci_j2_mott_analysis', \
            str(imp)]
    sp.call(cmd)

from pyglib.mbody.multiplets_analysis_soc import calc_save_atomic_states, \
        plot_atomic_states

# Calculate local histograms for impurity 1
# Check the 100 dominant eigen-states.
calc_save_atomic_states(imp=imp, num_ev=100)

# Plot the histogram
plot_atomic_states(imp=imp, num_label=5)
