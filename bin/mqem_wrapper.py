#!/usr/bin/env python
import numpy as np
import sys, subprocess,json,argparse, os,shutil
from scipy import interpolate

def gaux_sig(omega,sig):
    amat=np.zeros((4,4), dtype=complex)
    amat[0,:]=[1.0, 1.0/omega[-1], 1.0/omega[-1]**2, 1.0/omega[-1]**3]
    amat[1,:]=[1.0, -1.0/omega[-1], 1.0/omega[-1]**2, -1.0/omega[-1]**3]
    amat[2,:]=[1.0, 1.0/omega[-2], 1.0/omega[-2]**2, 1.0/omega[-2]**3]
    amat[3,:]=[1.0, -1.0/omega[-2], 1.0/omega[-2]**2, -1.0/omega[-2]**3]


    bmat=np.zeros(4, dtype=complex)
    bmat[0]=sig[-1]
    bmat[1]=np.conj(sig[-1])
    bmat[2]=sig[-2]
    bmat[3]=np.conj(sig[-2])

    coeff=np.linalg.solve(amat, bmat)
    # print('coeff', coeff)
    # gaux=(sig-np.real(coeff[0]))/np.real(coeff[1])
    
    
    # return gaux, np.real(coeff[0]),np.real(coeff[1])

    # gaux=sig-np.real(coeff[0])

    return coeff




def main(option):
    sigfile=option['sigfile']
    gaux_mode=option['gaux_mode']
    default_model=option['default_model']

    # emin_outer=option['emin_outer']
    # emax_outer=option['emax_outer']
    # emin_inner=option['emin_inner']
    # emax_inner=option['emax_inner']
    # egrid=option['egrid']
    blur=float(option['smearing'])
    interpolation=option['interpolation']    

    
    sigin=np.loadtxt(sigfile)
    sig=sigin[:,1::2]+sigin[:,2::2]*1j
    omega=sigin[:,0]
    nomega=np.shape(sig)[0]
    norb=np.shape(sig)[1]
    print(norb,nomega)
    print(sigfile)
    
    print('gaux_mode: ',gaux_mode)

    
    gaux=np.zeros(nomega, dtype=complex)
    if os.path.exists('realFreq_Sw.dat_1_1'):
        os.remove('realFreq_Sw.dat_1_1')

    emax=0
    for iorb in range(norb):
    # for iorb in range(1):    
        # print('--------------------------')
        # print('orb:  ' +  str(iorb+1))
        # print('--------------------------')

        tail=np.real(gaux_sig(omega*1j, sig[:,iorb]))
        a_center=tail[2]/tail[1]        
        a_width=np.sqrt(tail[3]/tail[1]-(tail[2]/tail[1])**2)
        if (gaux_mode =='sigc'):
            print('iorb '+str(iorb+1)+' center:', a_center, 'width:', a_width)
            emax=max(round(abs(a_center)+a_width*30), emax)            

        if (gaux_mode =='g'):
            tail=np.real(gaux_sig(omega*1j, 1.0/(omega*1j-(sig[:,iorb]-tail[0]))))
            a_center=tail[2]/tail[1]        
            a_width=np.sqrt(tail[3]/tail[1]-(tail[2]/tail[1])**2)
            print('iorb '+str(iorb+1)+' center:', a_center, 'width:', a_width)        
            emax=max(round(abs(a_center)+a_width*30), emax)
    print('emax:', emax)
    print('\n')        

        
    for iorb in range(norb):
    # for iorb in range(1):    
        print('--------------------------')
        print('orb:  ' +  str(iorb+1))
        print('--------------------------')
        tail=np.real(gaux_sig(omega*1j, sig[:,iorb]))
        if (gaux_mode=='sigc'):
            gaux=sig[:,iorb]
        elif(gaux_mode=='g'):
            gaux=1.0/(omega*1j-(sig[:,iorb]-tail[0]))
        f=open('gaux.dat', 'w')
        h=open('original_'+str(iorb+1), 'w')        
        for iomega in range(nomega):
            f.write('%5s %3s %3s %20.10f %20.10f \n' %(iomega, 0, 0, np.real(gaux[iomega]), np.imag(gaux[iomega])))
            h.write('%20.10f %20.10f %20.10f \n' %(omega[iomega], np.real(gaux[iomega]), np.imag(gaux[iomega])))            
        f.close()
        h.close()
        f=open('mqem.input.toml', 'w')
        # f.write('auxiliary_inverse_temp_range = [0.001, 500.0]\n')        
        f.write("inputFile = \"gaux.dat\"\n")
        f.write("NumOrbit = 1\n")
        f.write("inverse_temp =  "+str(np.pi/omega[0])+'\n')
        # f.write("Egrid = "+egrid+"\n")
        # f.write("EwinOuterRight = "+emax_outer+"\n")
        # f.write("EwinOuterLeft = "+emin_outer+"\n")
        # f.write("EwinInnerRight = "+emax_inner+"\n")
        # f.write("EwinInnerLeft = "+emin_inner+"\n")
        f.write("Egrid = 400\n")
        f.write("EwinOuterRight = "+str(emax)+"\n")
        f.write("EwinOuterLeft = -"+str(emax)+"\n")
        f.write("EwinInnerRight = 2\n")
        f.write("EwinInnerLeft = -2\n")
        
        # f.write("NumIter = 20\n")
        f.write("blur = "+str(blur)+"\n")        
        
        f.write("default_model=\""+default_model+"\"\n")      
        f.close()
        subprocess.call('julia  $MQEM/src/mem.jl', shell=True)
        if os.path.exists('realFreq_Sw.dat_1_1'):        
            gaux_real=np.loadtxt('realFreq_Sw.dat_1_1')
        else:
            print('maxent for orbital '+str(iorb+1)+' failed')
            sys.exit()
        freal_temp=gaux_real[:,0]
        ecutoff=min(15, emax/4)
        n1=np.argmin(abs(freal_temp+ecutoff))
        n2=np.argmin(abs(freal_temp-ecutoff))
        nf=n2-n1
        
        if (iorb ==0):
            sig_real=np.zeros((nf, 2*norb+1))
        if (gaux_mode=='sigc'):
            sig_real[:,0]=gaux_real[n1:n2,0]
            sig_real[:,2*iorb+1]=gaux_real[n1:n2,1]
            sig_real[:,2*iorb+2]=gaux_real[n1:n2,2]
        elif (gaux_mode=='g'):        
            sig_real[:,0]=gaux_real[n1:n2,0]
            sigout=gaux_real[n1:n2,0]-1.0/(gaux_real[n1:n2,1]+gaux_real[n1:n2,2]*1j)
            sig_real[:,2*iorb+1]=np.real(sigout)+tail[0]
            sig_real[:,2*iorb+2]=np.imag(sigout)
            # sig_real[:,2*iorb+1]=real(sigout)
            # sig_real[:,2*iorb+2]=imag(sigout)

            
        os.remove('realFreq_Sw.dat_1_1')
        os.remove('Sw_SOLVER.full_fromRetardedSw.dat_0_0')        
        os.remove('gaux.dat')        
        shutil.move('spectral_function_0_0.dat_model', 'gaux_spectra_model_'+str(iorb+1))
        shutil.move('spectral_function_0_0.dat', 'gaux_spectra_'+str(iorb+1))        
        shutil.move('information.out', 'info_'+str(iorb+1))
        shutil.move('reproduce_0_0.out', 'reproduced_'+str(iorb+1))
        print('\n')
        print('\n')
        print('\n')
        print('\n')
        print('\n')

    if (interpolation):
        np.savetxt('sig_realaxis_0.dat', sig_real, header='  ')        
        sig_real_bare=sig_real
        deltae=min(blur, 0.01)
        nn=int(min(ecutoff, 12)/deltae)
        sig_real=np.zeros((2*nn+1, 2*norb+1))
        xnew=(np.arange(2*nn+1)-nn)*deltae
        sig_real[:,0]=xnew
        for ii in range(2*norb):
            x=sig_real_bare[:,0]
            y=sig_real_bare[:,ii+1]
            f = interpolate.interp1d(x, y, kind='cubic')
            sig_real[:,ii+1]=f(xnew)
        np.savetxt('sig_realaxis.dat', sig_real, header='  ')                    
    else:
        np.savetxt('sig_realaxis.dat', sig_real, header='  ')
        
    

    print('Analytical continuation has finished:)')
    print('self-energ on real axis is in \"sig_realaxis.dat.\"')

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='call MQEM and return self-energy on real axis. The name of the output file will be sig_realaxis.dat')
    parser.add_argument('sigfile', action='store', help='self-energy file on imaginary axis')
    parser.add_argument('gaux_mode', action='store', nargs='?', default='g', help='The way to construct auxuliary green\'s function. sigc: \Sigma-\Sigma(omega=\inf), g: 1/(i\omega-Sigma).  default: g')
    parser.add_argument('default_model', action='store', nargs='?', default='g_mat', help='The way to construct auxuliary green\'s function. f:flat, g:gaussian, g_mat:matrix version of gaussian.  default: g_mat')
    parser.add_argument('smearing', action='store', nargs='?', default='0.01', help='smaring width.  default: 0.01')
    parser.add_argument('interpolation', action='store', nargs='?', default=True, help='smaring width.  default: True')    
    # parser.add_argument('emin_outer', action='store', nargs='?', default='-400', help='real frequency parameter to define the outer coarse energy grid.  default: -400')
    # parser.add_argument('emax_outer', action='store', nargs='?', default='400', help='real frequency parameter to define the outer coarse energy grid.  default: 400')
    # parser.add_argument('emin_inner', action='store', nargs='?', default='-0.001', help='real frequency parameter to define the inner dense energy grid.  default: -0.001')
    # parser.add_argument('emax_inner', action='store', nargs='?', default='0.001', help='real frequency parameter to define the inner dense energy grid.  default: 0.001')
    # parser.add_argument('egrid', action='store', nargs='?', default='2', help='real frequency parameter to define the inner dense energy grid.  default: 100')    

    option = parser.parse_args()
    # print(type(option))
    # print(vars(option))
    with open("mqem_wrapper_option.json", "w") as outfile:  
        json.dump(vars(option), outfile, indent=4) 

    main(vars(option))
