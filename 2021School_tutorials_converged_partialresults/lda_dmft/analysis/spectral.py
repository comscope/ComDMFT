import sys, os
from pylab import *
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np

dat0=reshape(loadtxt('../band/wannier_band_non_interpolated.dat'), [650, 9, 3], order='F')
dat1=reshape(loadtxt('../band/wannier_band_qp_interpolated.dat'), [650, 9, 3], order='F')
dat2=reshape(loadtxt('../realaxis/spectral.dat'), [650, 2401, 3], order='F')

close(1);figure(1);
pcolor(dat2[:,:,0], dat2[:,:,1], dat2[:,:,2], cmap='gray_r');colorbar();clim([0, 1])
plot(dat0[:,:,0], dat0[:,:,2], 'r-', dat1[:,:,0], dat1[:,:,2], 'b-');ylim([-5,1]);grid('on')
xticks(array([1, 119,202,285,387,489,591,650]), ('G','H','N','G','P','H','P','N'))
xlim([1, 650])
show()
