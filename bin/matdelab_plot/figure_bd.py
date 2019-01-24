"""
A class to draw the combined band structure - density of states plot

This class manages the interactions with the underlying band structure and
density of states plots.
"""

import matplotlib.pyplot   as plt
import matplotlib.gridspec as gridspec
from band_structure    import band
from density_of_states import dos

class figure_bd:

    def __init__(self,plotType,properties):
        """
        Create a figure for the given plotType

        If the plotType is 'band' then just generate a figure for the
        band structure plot.

        If the plotType is 'dos' then just generate a figure for the
        density of states plot.

        If the plotType is 'bd' then just generate a figure for both
        the band structure and the density of states plots.

        Other methods can be used to modify the parameters of these
        plots.

        Properties is a dictionary of properties that may be set for
        the plot.
        """
        self.fig = plt.figure(1)
        #plt.ion()
        if plotType == "band":
            self.fig.band = band(self.fig.add_subplot(111))
            self.fig.dos  = None
        elif plotType == "dos":
            self.fig.band = None
            self.fig.dos  = dos(self.fig.add_subplot(111))
        elif plotType == "bd":
            gs = gridspec.GridSpec(1,2,width_ratios=[4,1])
            if "width_ratios" in properties:
                gs.set_width_ratios(properties["width_ratios"])
            self.fig.band = band(self.fig.add_subplot(gs[0]))
            self.fig.dos  = dos(self.fig.add_subplot(gs[1]))

    def set_energy_range(self,emin,emax):
        """
        Set the energy ranges for the plots
        """
        if self.fig.band:
            self.fig.band.set_energy_range(emin,emax)
        if self.fig.dos:
            self.fig.dos.set_energy_range(emin,emax)

    def set_path(self,path):
        """
        Set the path along the symmetry points in the band structure
        """
        if self.fig.band:
            self.fig.band.set_path(path)

    def set_number_range(self,nmin,nmax):
        """
        Set the number of state ranges for the plots
        """
        if self.fig.dos:
            self.fig.dos.set_number_range(nmin,nmax)

    def add_band(self,xlist,ylist):
        """
        Add a band to the band structure plot
        """
        if self.fig.band:
            self.fig.band.add_band(xlist,ylist)

    def add_dos(self,xlist,ylist):
        """
        Add a DOS to the Density of States plot
        """
        if self.fig.dos:
            self.fig.dos.add_dos(xlist,ylist)
