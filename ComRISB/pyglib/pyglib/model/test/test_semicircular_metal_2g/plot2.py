import numpy as np

data_u1 = np.loadtxt('mu_u3.5.dat').T
data_u2 = np.loadtxt('mu_u5.0.dat').T


import matplotlib.pyplot as plt
plt.figure(1)
plt.plot(data_u1[0], data_u1[1], label='u=3.5')
plt.plot(data_u2[0], data_u2[1], label='u=5.0')
plt.legend()
plt.ylabel('$E_{tot}$')
plt.xlabel('$\mu$')
plt.figure(2)
plt.plot(data_u1[0], data_u1[2], label='u=3.5')
plt.plot(data_u2[0], data_u2[2], label='u=5.0')
plt.legend()
plt.ylabel('$n_f$')
plt.xlabel('$\mu$')
plt.show()
