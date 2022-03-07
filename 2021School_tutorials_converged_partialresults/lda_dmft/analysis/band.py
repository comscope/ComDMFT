import sys, os
from pylab import *
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np



dat0=reshape(loadtxt('../band/wannier_band_non_interpolated.dat'), [650, 9, 3], order='F')
dat1=reshape(loadtxt('../band/wannier_band_qp_interpolated.dat'), [650, 9, 3], order='F')
close(1);figure(1);plot(dat0[:,:,0], dat0[:,:,2], 'r-', dat1[:,:,0], dat1[:,:,2], 'b-');ylim([-5,1]);grid('on')
xticks(array([1, 119,202,285,387,489,591,650]), ('G','H','N','G','P','H','P','N'))
# plot(dat1[:,0]*30,dat1[:,1],'r+')
xlim([1, 650])
show()
# savefig('spectra_lda_dmft.png')
