# DMFT_MatDeLab bin/matdelab_plot

This directory contains a number of Python based components to plot
results from the various MatDeLab codes. At the moment the plots are 
related to the bandstructure and density of states (DOS) of GW calculations.

For the purpose of actually generating plots the scripts use the Matplotlib
library in Python (see http://matplotlib.org/).

## Band structure

The band structure is a representation of the one-electron state energies as
functions of k. It is customary to plot the band structure along a path 
following high-symmetry directions. The band structure plot allows one to 
establish whether the band gap of a material is a direct (HOMO-LUMO separation
at one k-vector) or indirect band gap (the lowest LUMO energy minus the highest
HOMO energy throughout the Brillouin zone).

## Density of States

The density of states gives a representation of the number of states as a
function of the energy throughout the Brillouin zone. Just considering
special directions is not sufficient, instead the whole Brillouin zone has to
be sampled evenly. Because the density of states involves the whole Brillouin
zone the information in this plot is in some sense more complete than that
of the band structure. 

## Common plots

There are three different common plots for displaying the information
discussed here. In all three different plot the Fermi level is chosen to
be at 0 energy. The energy is always plotted in the vertical direction.
Otherwise the plots are one of the three types below:

- the band structure 
- the density of states
- the band structure and the density of states

In the bandstructure plots the x-axis gives point along the high symmetry
paths. In the density of states the x-axis gives a count (normalized with the
Brillouin zone volume?) of the states at each energy level.
