"""
A class to draw a band structure plot.

This class uses Matplotlib underneath to implement most the functionality.
"""

import matplotlib.pyplot as plt
import array

class band:

    def __init__(self,axes):
        """
        Create a band structure plot by taking an axes instance and 
        modifying it in a suitable way. The axes instance is generated
        by a subplots command within a figure instance.
        """
        self.band = axes
        self.band.set_xlim(auto=True)
        self.band.axhline(y=0.0,color="black")

    def set_energy_range(self,emin,emax):
        """
        Set the energy ranges for the band structure.
        """
        self.band.set_ylim(bottom=emin,top=emax)

    def set_path(self,path):
        """
        Set the path along the symmetry points
        """
        self.band.path = path # store the path as that might come in handy
        (junk,xmin) = path[0]
        (junk,xmax) = path[-1]
        self.band.set_xlim(left=xmin,right=xmax)
        labels = array.array('c')
        xpos   = array.array('d')
        for xtic in path:
            (label,xcoord) = xtic
            labels.append(label)
            xpos.append(xcoord)
            self.band.axvline(x=xcoord,color="black")
        self.band.set_xticks(xpos,minor=False)
        self.band.set_xticklabels(labels,minor=False)

    def plot_band(self,*args,**kwargs):
        """
        Add and plot band structure data. The kwargs specify which data
        is to be plotted.
        """
        self.band.plot(kwargs)

    def add_band(self,xlist,ylist):
        """
        Add and plot a band to the band structure.
        """
        self.band.plot(xlist,ylist)
