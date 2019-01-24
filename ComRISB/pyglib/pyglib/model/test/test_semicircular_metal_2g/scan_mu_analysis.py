import numpy as np
from pyglib.model import circauxi
import shutil,subprocess
import h5py

mu_list = [1.569]+ \
        np.arange(1.57,1.6,0.01).tolist() + np.arange(1.6,3.0,0.1).tolist()
cmd = ['/home/ykent/WIEN_GUTZ/bin2/CyGutz', '-r', '-1']
e_list = []
n_list = []
u = 5.0
frec = open('mu_u5.0.dat', 'w')
for mu in mu_list:
    circauxi.gutz_model_setup(u=u, nmesh=5000, norb=3, tiny=0.0, mu=mu)
    shutil.copyfile('U{0}/WH_RL_INIT.h5_u{0}mu{1}'.format(u,mu),'WH_RL_INIT.h5')
    subprocess.call(cmd)
    shutil.copyfile('WH_RL_BEST.h5', 'WH_RL_INIT.h5_u{}mu{}'.format(u,mu))

    with h5py.File('EMBED_HAMIL_RES_1.h5', 'r') as f:
        n = f['/DM'][0,0] + f['/DM'][1,1]
        n_list.append(n)
    with h5py.File('GLOG.h5', 'r') as f:
        e = f['/etot_model'][0] + u/2*n
        e_list.append(e)

    frec.write('{} {} {}\n'.format(mu, e.real, n.real))
    frec.flush()
frec.close()
