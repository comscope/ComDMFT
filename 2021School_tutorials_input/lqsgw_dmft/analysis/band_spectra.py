import sys, os
from pylab import *
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np

dat1=reshape(loadtxt('../realaxis/wannier_band_non_interpolated.dat'), [84, 9, 3], order='F')
dat2=reshape(loadtxt('../realaxis/spectral.dat'), [84, 2401, 3], order='F')

close(1);figure(1);
pcolor(dat2[:,:,0], dat2[:,:,1], dat2[:,:,2], cmap='gray_r');colorbar();clim([0, 10])
plot(dat1[:,:,0], dat1[:,:,2], 'r-');ylim([-5,1]);grid('on')
xticks(array([1, 84]), ('N','G'))
xlim([1, 84])
show()
