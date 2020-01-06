#!/usr/bin/env python
from pyglib.iface.ifwannier import wannier_den_matrix
import sys


if "-w" in sys.argv:
    wannier_path = sys.argv[sys.argv.index("-w")+1]
else:
    wannier_path = "../wannier/"

wannier_den_matrix(wannier_path=wannier_path)
