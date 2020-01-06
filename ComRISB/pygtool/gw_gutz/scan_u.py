import numpy as np
import os
import subprocess
import shutil

j = 0.9
u_list = np.arange(2.0,7.0,1.0)

cmd_cygutz = ['/home/ykent/WIEN_GUTZ/bin2/CyGutz']

for u in u_list:
    print(' Running u = {}'.format(u))
    init_input = 'n\ny\nn\n1\n1\ny\ny\nd\n{} 0.9\nn\n-1\n-1'.format(u)
    with open('tmp.inp', 'w') as f:
        f.write(init_input)
    os.system('/home/ykent/WIEN_GUTZ/bin2/init_ga.py -m -u ev < "tmp.inp"')
    subprocess.call(cmd_cygutz)
    path = 'U{}J{}'.format(u,j)
    if not os.path.exists(path):
        os.mkdir(path)
    shutil.copy('GLOG.h5', path)
    shutil.copy('GBANDS_0.h5', path)
    shutil.copyfile('WH_RL_OUT.h5', path+'/WH_RL_INIT.h5')
    shutil.copyfile('WH_RL_OUT.h5', 'WH_RL_INIT.h5')
