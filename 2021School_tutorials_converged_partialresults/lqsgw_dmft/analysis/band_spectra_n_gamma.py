import sys, os
from pylab import *
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np

dat0=reshape(loadtxt('../band_n_gamma/wannier_band_non_interpolated.dat'), [202, 9, 3], order='F')
dat1=reshape(loadtxt('../band_n_gamma/wannier_band_qp_interpolated.dat'), [202, 9, 3], order='F')
dat2=reshape(loadtxt('../realaxis_n_gamma/spectral.dat'), [202, 2401, 3], order='F')

close(1);figure(1);
pcolor(dat2[:,:,0], dat2[:,:,1], dat2[:,:,2], cmap='gray_r');colorbar();clim([0, 10])
plot(dat0[:,:,0], dat0[:,:,2], 'r-', dat1[:,:,0], dat1[:,:,2], 'b-');ylim([-5,1]);grid('on')
xticks(array([1, 84,202]), ('N', 'G','H'))
xlim([1, 202])
show()
