import h5py
import numpy as np
from pyglib.model import circauxi
import shutil,subprocess

mu_list = np.arange(1.568,1.56,-10.001)
cmd = ['/home/ykent/WIEN_GUTZ/bin2/CyGutz', '-r', '-1']
u=5.0
flog = open('log', 'w')
for mu in mu_list:
    circauxi.gutz_model_setup(u=u, nmesh=5000, norb=3, tiny=0.0, mu=mu)
    subprocess.call(cmd)
    shutil.copyfile('WH_RL_BEST.h5', 'WH_RL_INIT.h5')
    shutil.copyfile('WH_RL_BEST.h5', 'WH_RL_INIT.h5_u{}mu{}'.format(u, mu))
    with h5py.File('GLOG.h5', 'r') as f:
        flog.write(' {}  {}\n'.format(mu, f['/rl_maxerr'][0]))
        flog.flush()

