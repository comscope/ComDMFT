import numpy
from pyglib.model.special import semicircular
import matplotlib.pyplot as plt


# Uniform e-mesh with different weight.
nmesh = 500
e_list = numpy.linspace(-1,1,nmesh)
sc = semicircular()
plt.plot(e_list, sc.dos(e_list), label='dos')
plt.plot(e_list, sc.cdos(e_list), label='cdos')
plt.xlim(-1.0, 1.0)
plt.ylim(0, 1.0)
plt.legend()
plt.xlabel('$\omega$')
plt.savefig('semicir_dos.png')
