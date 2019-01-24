#!/usr/bin/env python

import argparse
import numpy as np
import os, sys, shutil, subprocess, glob
import os.path
from numpy import pi
from scipy import *


def main(options):
    print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    print 
    print 'copy files for the real-axis calculation'
    print    

    dir_imag=os.path.abspath(options.lowh_directory)
    wan_dir=os.path.abspath(options.wan_directory)
    curdir=os.getcwd()
    shutil.copy(dir_imag+"/dc.dat", curdir)
    shutil.copy(dir_imag+"/ef.dat", curdir)
    shutil.copy(dir_imag+"/trans_basis.dat", curdir)
    shutil.copy(options.self_energy, curdir+'/sig.dat')
    files = glob.iglob(dir_imag+"/zinv_m1.dat")
    for filename in files:    
	shutil.copy(filename, curdir)
    os.system("cp "+dir_imag+"/wannier.dat"+" "+curdir)
    os.system("cp "+wan_dir+"/wannier.inip"+" "+curdir)
    if (options.mode == 3):
	print 'please generate kpath.dat file to calculate spectral function'


    print 'you should create job submission script'
	# print 	

    # if (options.ini_file):
    nomega=shape(loadtxt('sig.dat'))[0]
    print 'comlowh.ini file has been created'
    print 
    f=open(dir_imag+'/comlowh.ini', 'r')
    cnt=0
    numcixline=0			
    for cnt in range(5):
	line=f.readline()
	# 	    print line
	if (cnt ==2):
	    norb=map(int, line.split('#')[0].split())
	if (cnt ==3):
	    equval=map(int, line.split('#')[0].split())
	    maxequval=amax(equval)
	    minequval=amin(equval)		
	    for ii in range(minequval, maxequval+1):
		for jj in range(shape(norb)[0]):
		    if equval[jj] ==ii:
			numcixline=numcixline+norb[jj]
			break
    f.close()
	
    f=open(dir_imag+'/comlowh.ini', 'r')
    inifile=f.readlines()
    inifile[0]=str(options.mode)+'\n'# +' #'+inifile[0].split('#')[1]
    inifile[6+numcixline]=str(nomega)+'\n'# +' #'+inifile[8+numcixline].split('#')[1]
    inifile[8+numcixline]=str(options.broadening)+'\n'# +' #'+inifile[10+numcixline].split('#')[1]
    if ((options.kmesh_b1_for_dos!=None) & (options.kmesh_b2_for_dos!=None) & (options.kmesh_b3_for_dos!=None)):
	inifile[14+numcixline]=options.kmesh_b1_for_dos+'  '+options.kmesh_b2_for_dos+'  '+options.kmesh_b3_for_dos+'\n'
    f.close()
    # 	print inifile
    g=open('comlowh_new.ini', 'w')
    g.writelines(inifile)
    g.close()
    shutil.move("comlowh_new.ini", "comlowh.ini")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='prepare inputs of comlowh calculation on real axis')
    parser.add_argument('broadening', action='store', help='broadening')        
    parser.add_argument('lowh_directory', action='store', help='lowh directory')
    parser.add_argument('wan_directory', action='store', help='wannier directory')
    parser.add_argument('self_energy', action='store', help='real-axis self-energy')
    parser.add_argument('-m', '--mode', action='store', type=int, default=3, help='If 3, code calculates spectral function along the high symmetry line defined in \'kpath.dat\'. If it is 2, it calculates projected density of states. Default: 3')
    parser.add_argument('kmesh_b1_for_dos', action='store', nargs='?', help='finer kmesh along b1 axis for the DOS. Optional')
    parser.add_argument('kmesh_b2_for_dos', action='store', nargs='?', help='finer kmesh along b2 axis for the DOS. Optional')
    parser.add_argument('kmesh_b3_for_dos', action='store', nargs='?', help='finer kmesh along b3 axis for the DOS. Optional')    
    
    
    
    options = parser.parse_args()
#     print options.wandirectory
    
    main(options)
