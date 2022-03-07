import sys, os
from pylab import *
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np



delta=loadtxt('../delta.dat')
close(1);figure(1);plot(delta[:,0],delta[:,1], 'r');xlim([0, 20]);grid('on')
xlabel('omega(eV)')
ylabel('Re delta(eV)')

close(2);figure(2);plot(delta[:,0],delta[:,2], 'r');xlim([0, 20]);grid('on')

xlabel('omega(eV)')
ylabel('Im delta(eV)')
show()
