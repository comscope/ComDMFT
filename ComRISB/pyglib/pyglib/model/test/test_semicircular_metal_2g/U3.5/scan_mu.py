import numpy as np
from pyglib.model import circauxi
import shutil,subprocess

mu_list = np.arange(2.1,3.0,10.1)
cmd = ['/home/ykent/WIEN_GUTZ/bin2/CyGutz', '-r', '-1']
u=3.5
for mu in mu_list:
    circauxi.gutz_model_setup(u=u, nmesh=5000, norb=3, tiny=0.0, mu=mu)
    subprocess.call(cmd)
    shutil.copyfile('WH_RL_BEST.h5', 'WH_RL_INIT.h5')
    shutil.copyfile('WH_RL_BEST.h5', 'WH_RL_INIT.h5_u{}mu{}'.format(u, mu))

