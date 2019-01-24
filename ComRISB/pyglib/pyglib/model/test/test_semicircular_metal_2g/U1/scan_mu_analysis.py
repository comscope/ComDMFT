import numpy as np
from pyglib.model import circauxi
import shutil,subprocess
import h5py

mu_list = np.arange(0.0,1.0,0.1)
cmd = ['/home/ykent/WIEN_GUTZ/bin2/CyGutz', '-r', '-1']
e_list = []
n_list = []
u = 1.0
frec = open('mu.dat', 'w')
for mu in mu_list:
    circauxi.gutz_model_setup(u=1.0, nmesh=5000, norb=3, tiny=0.0, mu=mu)
    shutil.copyfile('WH_RL_INIT.h5_u1mu{}'.format(mu),'WH_RL_INIT.h5')
    subprocess.call(cmd)
    shutil.copyfile('WH_RL_BEST.h5', 'WH_RL_INIT.h5_u1mu{}'.format(mu))

    with h5py.File('EMBED_HAMIL_RES_1.h5', 'r') as f:
        n = f['/DM'][0,0] + f['/DM'][1,1]
        n_list.append(n)
    with h5py.File('GLOG.h5', 'r') as f:
        e = f['/etot_model'][0] + u/2*n
        e_list.append(e)

    frec.write('{} {} {}\n'.format(mu, e.real, n.real))
    frec.flush()
frec.close()
