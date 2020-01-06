#!/usr/bin/env python

import sys, glob
from pyglib.gutz.init import initialize
from pyglib.iface.ifwien import h4set_indmfl


initialize()
fstruct = glob.glob('*struct')
if '-no-indmfl' not in sys.argv and len(fstruct) == 1:
    h4set_indmfl()
