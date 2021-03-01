#!/usr/bin/env python

import argparse
import sys, os, shutil
import numpy as np
import subprocess
from pylab import *

flapw= loadtxt(sys.argv[1])
wannier= loadtxt(sys.argv[2])

plot(flapw[:,0]/amax(flapw[:,0]), flapw[:,1], 'r.', wannier[:,0]/amax(wannier[:,0]), wannier[:,4], 'b.')



