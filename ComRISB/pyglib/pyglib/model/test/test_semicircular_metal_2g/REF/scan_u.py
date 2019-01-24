import numpy as np
from pyglib.model import circauxi
import shutil,subprocess


cmd = ['/home/ykent/WIEN_GUTZ/bin2/CyGutz', '-r', '-1']
for i,u in enumerate(np.arange(1.0,0.9,-10)):
    print(' Running with u = {}'.format(u))
    circauxi.gutz_model_setup(u=u, nmesh=5000, norb=3, tiny=0.0, mu=0.0)
    subprocess.call(cmd)
    shutil.copyfile('WH_RL_BEST.h5', 'WH_RL_INIT.h5')
    shutil.copyfile('WH_RL_BEST.h5', 'WH_RL_INIT.h5_{}'.format(u))

