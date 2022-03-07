import sys, os
from pylab import *
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np



sig=loadtxt('../sig.dat')
close(1);figure(1);plot(sig[:,0],sig[:,1], 'r');xlim([0, 20]);grid('on')
xlabel('omega(eV)')
ylabel('Re Sigma(eV)')
close(2);figure(2);plot(sig[:,0],sig[:,2], 'r');xlim([0, 20]);grid('on')
xlabel('omega(eV)')
ylabel('Im Sigma(eV)')
show()
