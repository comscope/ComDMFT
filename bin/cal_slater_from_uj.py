#!/usr/bin/env python
import argparse
import numpy as np
import os, sys, shutil, subprocess, glob
import os.path
from numpy import pi
from scipy import *

def main(options):

    problem=options.subshell
    uval=float(options.uval)
    jval=float(options.jval)    

    if (problem=='s'):
	F0=uval
	print 'F0=', F0

    elif (problem=='p'):
	F0=uval
	F2=5*jval
	print 'F0=', F0
	print 'F2=', F2
	
    elif (problem=='d'):
	F0=uval
	F2=14.0/1.625*jval
	F4=8.75/1.625*jval
	print 'F0=', F0
	print 'F2=', F2	
	print 'F4=', F4

    elif (problem=='f'):
	F0=uval
	F2=6435.0/(286+195*0.668+250*0.494)*jval
	F4=0.668*F2
	F6=0.494*F2
	print 'F0=', F0
	print 'F2=', F2	
	print 'F4=', F4
	print 'F6=', F6		


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='calculate slater-condon parameters with inputs of U and J') 
    parser.add_argument('subshell', action='store', help='subshell either s, p, d or f')
    parser.add_argument('uval', action='store', help='U value in eV')
    parser.add_argument('jval', action='store', help='J value in eV')    
    options = parser.parse_args()
    
    main(options)	
