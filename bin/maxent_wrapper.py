#!/usr/bin/env python
# developed by Siheon Ryee (siheonryee@gmail.com) and Sangkook Choi(sangkookchoi@gmail.com)

import argparse
import sys, os, shutil
import numpy as np
import subprocess
pi = 3.1415926535897932386
kB = 0.00008617328149741


def write_params(temperature,error):
    Ntau=int(300/float(temperature)*3000)
    f=open('maxent_params.dat', 'w')	
    f.write("params={'statistics': 'fermi', # fermi/bose\n")
    f.write("        'Ntau'      : "+str(Ntau)+",     # Number of time points\n")
    f.write("        'L'         : 70.0,     # cutoff frequency on real axis\n")
    f.write("        'Nw'        : 501,     # number of frequency points on real axis\n")
    f.write("        'gwidth'    : 140.0,  # width of gaussian\n")
    f.write("        'idg'       : 1,       # error scheme: idg=1 -> sigma=deltag ; idg=0 -> sigma=deltag*G(tau)\n")
    f.write("        'deltag'    : "+str(error)+",   # error\n")
    f.write("        'Asteps'    : 4000,    # anealing steps\n")
    f.write("        'alpha0'    : 1000,    # starting alpha\n")
    f.write("        'x0'        : 0.01,    # low energy cutfoff\n")    
    f.write("        'min_ratio' : 0.001,    # condition to finish, what should be the ratio\n")
    f.write("        'iflat'     : 1,       # iflat=0 : constant model, iflat=1 : gaussian of width gwidth, iflat=2 : input using file model.dat\n")
    f.write("        'Nitt'      : 1000,    # maximum number of outside iterations\n")
    f.write("        'Nr'        : 0,       # number of smoothing runs\n")
    f.write("        'bwdth'     : 0.03,    # smoothing width\n")
    f.write("        'Nf'        : 5,      # to perform inverse Fourier, high frequency limit is computed from the last Nf points\n")
    f.write("         }\n")	
    return None


def main(option):
    #if len(sys.argv)<2:
    #  print 'give input file sig.dat'
    #  sys.exit(0)
    
    #sigfile = sys.argv[1]
    sigfile = option.sig

    sigdata = np.loadtxt(sigfile).transpose()
    nb2=0; nz2=0;
    for i in range(1,len(sigdata)):
      if sum(abs(sigdata[i]))>0:
        nb2 +=1
      else:
        nz2 +=1
    nb = nb2/2
    nz = nz2/2
    
    
    # make temporary self-energy 
    sigtemp=[]
    sigtemp.append(sigdata[0])
    for b in range(0,nb):
      siginf = np.full(len(sigdata[1+2*b]),sigdata[1+2*b][-1])
      sigtemp.append(sigdata[1+2*b]-siginf)
      sigtemp.append(sigdata[2+2*b])
    sigtemp = np.array(sigtemp)
    
    np.savetxt('sig.temp', sigtemp.transpose())
    line = "#\n#"
    with open('sig.temp', 'r+') as f:
        file_data = f.read()
        f.seek(0,0)
        f.write(line.rstrip('\r\n') + '\n' + file_data)
    f.close()
    
    beta = pi/(sigdata[0][0])
    temp = (1.0/(beta*kB))
    if not os.path.isfile("maxent_params.dat"):
       write_params(temp,option.error)
    
    
    # path to Haule's maxent_run.py
    maxentpath = os.environ.get('WIEN_DMFT_ROOT')
    # run maxent_run.py
    cmd = maxentpath+'/maxent_run.py sig.temp'
    subprocess.call(cmd, shell=True)
    
    
    # collect and re-organize
    sig_ac = np.loadtxt('Sig.out').transpose()
    sigout=[]
    sigout.append(sig_ac[0])
    for b in range(0,nb):
      siginf = np.full(len(sig_ac[0]),sigdata[1+2*b][-1])
      sigout.append(sig_ac[1+2*b]+siginf)
      sigout.append(sig_ac[2+2*b])
    sigout = np.array(sigout)
    
    # The final output
    np.savetxt('sig_realaxis.dat', sigout.transpose())
    line = "#"
    with open('sig_realaxis.dat', 'r+') as f:
        file_data = f.read()
        f.seek(0,0)
        f.write(line.rstrip('\r\n') + '\n' + file_data)
    f.close()
    os.remove('Sig.out')
    os.remove('sig.temp')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='call maxent_run.py from EDMFTF package and return self-energy on real axis. If maxent_params.dat is not present in the directory, it generates one. The name of the output file will be sig_realaxis.dat')
    parser.add_argument('sig', action='store', help='self-energy file on imaginary axis')
    parser.add_argument('error', action='store', nargs='?', default='0.05', help='Errors for the maxent.Optional. Default value=0.05' )   

    option = parser.parse_args()
    main(option)
