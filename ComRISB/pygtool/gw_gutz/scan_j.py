import numpy as np
import os
import subprocess
import shutil

u = 4.8
j_list = np.arange(0.0,1.1,0.1)

cmd_cygutz = ['/home/ykent/WIEN_GUTZ/bin2/CyGutz']

for j in j_list:
    print(' Running j = {}'.format(j))
    init_input = 'n\ny\nn\n1\n1\ny\ny\nd\n4.8 {}\nn\n-1\n-1'.format(j)
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
