import sys, os
from pylab import *
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np



sig=loadtxt('../sig.dat')
dc=loadtxt('../sig_dc.dat')
close(1);figure(1);plot(sig[:,0],sig[:,1], 'r', dc[:,0], dc[:,1], 'b-');xlim([0, 20]);grid('on')
xlabel('omega(eV)')
ylabel('Re Sigma(eV)')
legend(['imp', 'dc'])
close(2);figure(2);plot(sig[:,0],sig[:,2], 'r', dc[:,0], dc[:,2], 'b-');xlim([0, 20]);grid('on')
legend(['imp', 'dc'])
xlabel('omega(eV)')
ylabel('Im Sigma(eV)')
show()
