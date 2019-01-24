import numpy as np
from pyglib.model import circauxi
import shutil,subprocess

u_list = np.arange(5.0, 5.1, 0.1)
cmd = ['/home/ykent/WIEN_GUTZ/bin2/CyGutz', '-r', '-1']
mu=2.1
for u in u_list:
    circauxi.gutz_model_setup(u=u, nmesh=5000, norb=3, tiny=0.0, mu=mu)
    subprocess.call(cmd)
    shutil.copyfile('WH_RL_BEST.h5', 'WH_RL_INIT.h5')
    shutil.copyfile('WH_RL_BEST.h5', 'WH_RL_INIT.h5_u{}mu{}'.format(u, mu))

