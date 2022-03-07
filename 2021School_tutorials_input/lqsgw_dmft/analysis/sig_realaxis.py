import sys, os
from pylab import *
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np

dat0=loadtxt('/home/max/codes/Compiled_ComsuiteCode/ComsuiteV2/tutorials_converged/lqsgw_dmft/maxent/sig_realaxis.dat')

close(1);figure(1);
plot(dat0[:,0], dat0[:,1], 'r-');
xlim([-10, 10]);grid('on')
xlabel('omega(eV)')
ylabel('Re Sigma(eV)')
close(2);figure(2);
plot(dat0[:,0],dat0[:,2], 'r-');
xlim([-10, 10]);grid('on')
xlabel('omega(eV)')
ylabel('Im Sigma(eV)')
show()
