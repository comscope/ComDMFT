"""
A class to draw a density of states plot.

This class uses Matplotlib underneath to implement most the functionality.
"""

import matplotlib.pyplot as plt

class dos:

    def __init__(self,axes):
        """
        Create a density of states plot by taking an axes instance and 
        modifying it in a suitable way. The axes instance is generated
        by a subplots command within a figure instance.
        """
        self.dos = axes
        self.dos.set_xlim(auto=True)
        self.dos.axhline(y=0.0,color="black")

    def set_energy_range(self,emin,emax):
        """
        Set the energy ranges for the DOS.
        """
        self.dos.set_ylim(bottom=emin,top=emax)

    def set_number_range(self,nmin,nmax):
        """
        Set the number of state ranges for the DOS.
        """
        self.dos.set_xlim(left=nmin,right=nmax)

    def plot_dos(self,*args,**kwargs):
        """
        Add and plot density of states data. The kwargs specify which data
        is to be plotted.
        """
        self.dos.plot(kwargs)

    def add_dos(self,elist,dlist):
        """
        Add and plot density of states data. 
        - elist is the list of energies
        - dlist is the density of states at each energy
        """
        self.dos.plot(dlist,elist)
