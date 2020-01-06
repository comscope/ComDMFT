import numpy as np

data_u1 = np.loadtxt('mu_u1.dat').T
data_u25 = np.loadtxt('mu_u2.5.dat').T


import matplotlib.pyplot as plt
plt.figure(1)
plt.plot(data_u1[0], data_u1[1], label='u=1.0')
plt.plot(data_u25[0], data_u25[1], label='u=2.5')
plt.legend()
plt.ylabel('$e_{tot}$')
plt.xlabel('$\mu$')
plt.figure(2)
plt.plot(data_u1[0], data_u1[2], label='u=1.0')
plt.plot(data_u25[0], data_u25[2], label='u=2.5')
plt.legend()
plt.ylabel('$n_f$')
plt.xlabel('$\mu$')
plt.show()
