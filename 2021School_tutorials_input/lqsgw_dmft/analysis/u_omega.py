import sys, os
from pylab import *
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np


uu=loadtxt('../u_slater.dat')
ww=loadtxt('../w_slater.dat')
vv=loadtxt('../v_slater.dat')
close(1);figure(1);plot(uu[:,0],uu[:,1], 'r', ww[:,0],ww[:,1], 'b', [0, 3000], array([1,1])*vv, 'k--');xlim([0, 40]);grid('on')
ylabel('U(eV)')
xlabel('nu(eV)')
legend(['u', 'w', 'v'])
#savefig('delta_real.png')


show()
