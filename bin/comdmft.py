#!/usr/bin/env python
from __future__ import print_function
import argparse
import numpy as np
import os, sys, shutil, subprocess, glob
import re
from numpy import pi
from scipy import *
import json
from tabulate import tabulate
from itertools import chain
# from scipy.interpolate import interp1d
# trans_basis_mode: 0, use wannier function as basis set
# trans_basis_mode: 1, use transformation matrix to rotate the basis set. this matrix doesn't change as a function of iteration.
# trans_basis_mode: 2, use transformation matrix to rotate the basis set. this matrix does change as a function of iteration. this matrix diagonalize the spectral function at the chemical potential.


def open_h_log(control):
    control['h_log']=open('./cmd.log', 'a')
    print('',                                 file=control['h_log'],flush=True)
    print('*********************************',file=control['h_log'],flush=True)
    print('             ComDMFT',             file=control['h_log'],flush=True)
    print('*********************************',file=control['h_log'],flush=True)
    print('',                                 file=control['h_log'],flush=True)
    #DEBUG
    control['h_log'].flush()
    os.fsync(control['h_log'].fileno())
    #DEBUG


    return None

def close_h_log(control):
    control['h_log'].close()
    return None

def read_comdmft_ini():
    vglobl={}
    vlocal={}
    if sys.version_info.major<3:
        execfile('comdmft.ini',vglobl,vlocal)
    else: # sys.version_info.major==3
        with open('comdmft.ini') as f_ini:
            code = compile(f_ini.read(), "comdmft.ini", 'exec')
            exec(code, vglobl, vlocal)
            f_ini.close()
    # print vglobl
    # print 'here'
    control=vlocal['control']
    wan_hmat=vlocal['wan_hmat']
    imp=vlocal['imp']

    control['name']='control'
    wan_hmat['name']='wan_hmat'
    imp['name']='imp'



    open_h_log(control)        

    control['comsuitedir']=os.environ.get('COMSUITE_BIN')
    if not control['comsuitedir']:
        print("Error: Environment variable COMSUITE_BIN is not defined.", file=control['h_log'],flush=True)
        sys.exit()
    print('comsuitedir', control['comsuitedir'])
    control['conv_table']=[]

    ###    in control

    control['top_dir']=os.path.abspath('./')        
    check_key_in_string('method', control, control)
    control['sigma_mix_ratio']=control.get('sigma_mix_ratio', 0.5)

    control['doping']=control.get('doping', 0.0)
    
    control['u_mode']=control.get('u_mode', 'bnse')

    control['trans_basis_mode']=control.get('trans_basis_mode', 0)
    if (control['trans_basis_mode']==1):
        check_key_in_string('trans_basis', control,control)
    elif (control['trans_basis_mode']==2):
        check_key_in_string('metal_threshold', control,control)

        
        
    check_key_in_string('spin_orbit', control, control)
    check_key_in_string('impurity_problem', control, control)    
    check_key_in_string('impurity_problem_equivalence', control, control)
    check_key_in_string('initial_lattice_dir', control, control)
    control['initial_lattice_dir']=os.path.abspath(control['initial_lattice_dir'])        
    control['allfile']=find_allfile(control['initial_lattice_dir'])
    if ('dc_directory' not in control):
        control['dc_directory']='./dc'
    control['dc_directory']=os.path.abspath(control['dc_directory'])
    if ('impurity_directory' not in control):
        control['impurity_directory']='./impurity'
    control['impurity_directory']=os.path.abspath(control['impurity_directory'])
    if ('lowh_directory' not in control):
        control['lowh_directory']='./lowh'
    control['lowh_directory']=os.path.abspath(control['lowh_directory'])
    if ('wannier_directory' not in control):
        control['wannier_directory']='./wannier'
    control['wannier_directory']=os.path.abspath(control['wannier_directory'])    

    if ('initial_self_energy' in control):
        control['initial_self_energy'] =os.path.abspath(control['initial_self_energy'])
        if (control['trans_basis_mode']!=0):
            check_key_in_string('trans_basis', control,control)

            
    if ('dc_mat_to_read' in control):
        control['dc_mat_to_read'] =os.path.abspath(control['dc_mat_to_read'])

    if (control['method']=='lda+dmft'):
        control['convergence_header']=['step','i_outer','i_latt','i_imp','causality','delta_rho','w_sp_min','w_sp_max', 'mu', 'std_sig', 'n_imp', 'histo_1', 'histo_2', 'ctqmc_sign']
    if (control['method']=='lqsgw+dmft'):
        control['convergence_header']=['step','i_imp','causality','static_f0','w_sp_min','w_sp_max', 'mu', 'std_sig', 'n_imp', 'histo_1', 'histo_2', 'ctqmc_sign']

    # mpi_prefix
    if ('mpi_prefix' in control):
        control['mpi_prefix_lowh']=control.get('mpi_prefix_lowh', control['mpi_prefix'])
        control['mpi_prefix_impurity']=control.get('mpi_prefix_impurity', control['mpi_prefix'])
        control['mpi_prefix_wannier']=control.get('mpi_prefix_wannier', control['mpi_prefix'])
        if (control['method']=='lda+dmft'):
            control['mpi_prefix_lattice']=control.get('mpi_prefix_lattice', control['mpi_prefix'])
        if (control['method']=='lqsgw+dmft'):    
            control['mpi_prefix_dc']=control.get('mpi_prefix_dc', control['mpi_prefix'])

    # mpi_prefix_coulomb
    if ('mpi_prefix_coulomb' in control):
        check_key_in_string('nproc_k_coulomb', control, control)
        check_key_in_string('nproc_tau_coulomb', control, control)
    else:
        # temp=[int(x) for x in np.loadtxt(control['initial_lattice_dir']+'/k_tau_freq.dat')]
        temp=list(map(int,np.loadtxt(control['initial_lattice_dir']+'/k_tau_freq.dat')))
        control['mpi_prefix_coulomb'], control['nproc_k_coulomb'],control['nproc_tau_coulomb']=optimized_nproc_for_comcoulomb(control['mpi_prefix'], temp[0], temp[1],temp[2],temp[3])
        # print('mpi_prefix_coulomb', control['mpi_prefix_coulomb'], file=control['h_log'],flush=True)

    # max iteration
    if (control['method']=='lda+dmft'):
        control['max_iter_num_impurity']=control.get('max_iter_num_impurity', 1)
        control['max_iter_num_outer']=control.get('max_iter_num_outer', 50)
    elif (control['method']=='lqsgw+dmft'):
        control['max_iter_num_impurity']=control.get('max_iter_num_impurity', 50)


    # directory_name
    if (control['method']=='lda+dmft'):
        if ('lattice_directory' not in control):
            control['lattice_directory']='./lattice'
        control['lattice_directory']=os.path.abspath(control['lattice_directory'])
    if (control['method']=='lqsgw+dmft'):
        if  ('coulomb_directory' not in control):
            control['coulomb_directory']='./coulomb'
        control['coulomb_directory']=os.path.abspath(control['coulomb_directory'])


    if (control['method']=='lqsgw+dmft'):
        control['do_wannier']=True
        control['do_coulomb']=True
        control['do_dc']=True
        control['iter_num_impurity']=1
        control['iter_num_outer']=1
    elif (control['method']=='lda+dmft'):
        control['iter_num_outer']=1
        control['iter_num_impurity']=0

    control['restart']=control.get('restart', False)

    if (control['restart']):
        find_place_to_restart(control)
        if (control['method']=='lqsgw+dmft'):
            print('do_wannier', control['do_wannier'], file=control['h_log'],flush=True)
            print('do_coulomb', control['do_coulomb'], file=control['h_log'],flush=True)
            print('do_dc',      control['do_dc'],      file=control['h_log'],flush=True)

#   in wan_hmat
    check_key_in_string('kgrid', wan_hmat, control)
    check_key_in_string('froz_win_min', wan_hmat, control)
    check_key_in_string('froz_win_max', wan_hmat, control)
    wan_hmat['dis_win_min']=wan_hmat.get('dis_win_min', wan_hmat['froz_win_min'])
    wan_hmat['dis_win_max']=wan_hmat.get('dis_win_max', wan_hmat['froz_win_max']+40.0)
    control['proj_win_min']=control.get('proj_win_min', wan_hmat['dis_win_min'])
    control['proj_win_max']=control.get('proj_win_max', wan_hmat['dis_win_max'])

    wan_hmat['num_iter']=wan_hmat.get('num_iter', 0)
    wan_hmat['dis_num_iter']=wan_hmat.get('dis_num_iter', 100)


    wan_hmat['cut_low']=wan_hmat.get('cut_low', 0.4)
    wan_hmat['cut_froz']=wan_hmat.get('cut_froz', 0.10)
    wan_hmat['cut_total']=wan_hmat.get('cut_total', 0.0)        


    if (control['method']=='lqsgw+dmft'):
        wan_hmat['rmode']=wan_hmat.get('rmode', 0)
        wan_hmat['radfac']=wan_hmat.get('radfac', 1.0)        
    if (control['method']=='lda+dmft'):
        wan_hmat['rmode']=wan_hmat.get('rmode', 0)                
        wan_hmat['radfac']=wan_hmat.get('radfac', 1.0)        
    # in imp

    check_key_in_string('temperature', imp, control)    
    imp['beta']=1.0/(8.6173303*10**-5*imp['temperature'])

    if ('initial_self_energy' in control):
        control['n_omega']=np.shape(np.loadtxt(control['initial_self_energy']))[0]
    else:
        control['n_omega']=int(300.0/(2*pi/imp['beta']))

    control['omega']=(arange(control['n_omega'])*2+1)*pi/imp['beta']

    for key, value in imp.items():

        if (not (isinstance(imp[key], dict))):
            continue
        imp[key]['name']=key


        # imp[key]['para']=True
        # for ktemp in control['impurity_problem_equivalence'] :
        #     if (ktemp == -1):
        #         imp[key]['para']=False


        if (-1*int(key) in control['impurity_problem_equivalence']):
            imp[key]['para']=False
        else:
            imp[key]['para']=True

        imp[key]['problem']=control['impurity_problem'][control['impurity_problem_equivalence'].index(int(key))][1]
        if (control['method']=='lda+dmft'):
            check_key_in_string('f0', imp[key], control)
            check_key_in_string('f2', imp[key], control)
            check_key_in_string('f4', imp[key], control)
            if (imp[key]['problem']=='f'):
                check_key_in_string('f6', imp[key], control)    
        # elif (control['method']=='lqsgw+dmft'):
        #     check_key_in_string('boson_low_truncation', imp[key], control)

        check_key_in_string('thermalization_time', imp[key], control)
        check_key_in_string('measurement_time', imp[key], control)
        check_key_in_string('impurity_matrix', imp[key], control)
        if (control['trans_basis_mode']<2):
            imp[key]['impurity_matrix']=array(imp[key]['impurity_matrix'])
        else:
            print("impurity_matrix reset", file=control['h_log'],flush=True)
            nimp_orb=len(imp[key]['impurity_matrix'])
            imp[key]['impurity_matrix']=np.zeros((nimp_orb,nimp_orb), dtype='int')
            for ii in range(nimp_orb):
                imp[key]['impurity_matrix'][ii,ii]=ii+1

        print('here',                            file=control['h_log'],flush=True)
        print(type(imp[key]['impurity_matrix']), file=control['h_log'],flush=True)

        print(imp[key]['impurity_matrix'],       file=control['h_log'],flush=True)
        print('here',                            file=control['h_log'],flush=True)

        if (control['method']=='lda+dmft'):
            check_key_in_string('nominal_n', imp[key], control)

        check_key_in_string('green_cutoff', imp[key], control)
        imp[key]['susceptibility_cutoff']=imp[key].get('susceptibility_cutoff', 50)
        imp[key]['susceptibility_tail']=imp[key].get('susceptibility_tail', 300)        

        if ('coulomb' not in imp[key]):
            imp[key]["coulomb"]='full'



    control['sig_header']=['# omega(eV)']
    for ii in sorted(set(control['impurity_problem_equivalence'])):
        for jj in sorted(set(imp[str(abs(ii))]['impurity_matrix'].flatten().tolist())-{0}):
            control['sig_header'].append("Re Sig_{"+str(ii)+','+str(jj)+'}(eV)')
            control['sig_header'].append("Im Sig_{"+str(ii)+','+str(jj)+'}(eV)')    

    # print
    print('top_dir', control['top_dir'],                         file=control['h_log'],flush=True)
    if (control['method']=='lda+dmft'):
        print('lattice_directory', control['lattice_directory'], file=control['h_log'],flush=True)
    elif (control['method']=='lqsgw+dmft'):    
        print('coulomb_directory', control['coulomb_directory'], file=control['h_log'],flush=True)
    print('wannier_directory', control['wannier_directory'],     file=control['h_log'],flush=True)
    print('dc_directory', control['dc_directory'],               file=control['h_log'],flush=True)
    print('impurity_directory', control['impurity_directory'],   file=control['h_log'],flush=True)
    print('lowh_directory', control['lowh_directory'],           file=control['h_log'],flush=True)



    return control,wan_hmat,imp

def find_impurity_wan(control, wan_hmat):
    num_wann=np.shape(wan_hmat['basis'])[0]
    control['impurity_wan']=[]
    for ip in range(np.shape(control['impurity_problem'])[0]):
        if (control['spin_orbit']):
            if (control['impurity_problem'][ip][1].lower()=='f'):
                control['impurity_wan'].append([0]*14)
                for iwan in range(num_wann):
                    if ((wan_hmat['basis'][iwan]['atom']==control['impurity_problem'][ip][0]) and (wan_hmat['basis'][iwan]['l']==3)):
                        if (int(wan_hmat['basis'][iwan]['i']*2)==-1):
                            if (int(wan_hmat['basis'][iwan]['m']*2)==-5):
                                control['impurity_wan'][ip][0]=wan_hmat['basis'][iwan]['ind']
                            elif (int(wan_hmat['basis'][iwan]['m']*2)==-3):
                                control['impurity_wan'][ip][1]=wan_hmat['basis'][iwan]['ind']
                            elif (int(wan_hmat['basis'][iwan]['m']*2)==-1):
                                control['impurity_wan'][ip][2]=wan_hmat['basis'][iwan]['ind']
                            elif (int(wan_hmat['basis'][iwan]['m']*2)==1):
                                control['impurity_wan'][ip][3]=wan_hmat['basis'][iwan]['ind']
                            elif (int(wan_hmat['basis'][iwan]['m']*2)==3):
                                control['impurity_wan'][ip][4]=wan_hmat['basis'][iwan]['ind']
                            elif (int(wan_hmat['basis'][iwan]['m']*2)==5):
                                control['impurity_wan'][ip][5]=wan_hmat['basis'][iwan]['ind']
                        elif (int(wan_hmat['basis'][iwan]['i']*2)==1):
                            if (int(wan_hmat['basis'][iwan]['m']*2)==-7):
                                control['impurity_wan'][ip][6]=wan_hmat['basis'][iwan]['ind']
                            elif (int(wan_hmat['basis'][iwan]['m']*2)==-5):
                                control['impurity_wan'][ip][7]=wan_hmat['basis'][iwan]['ind']
                            elif (int(wan_hmat['basis'][iwan]['m']*2)==-3):
                                control['impurity_wan'][ip][8]=wan_hmat['basis'][iwan]['ind']
                            elif (int(wan_hmat['basis'][iwan]['m']*2)==-1):
                                control['impurity_wan'][ip][9]=wan_hmat['basis'][iwan]['ind']
                            elif (int(wan_hmat['basis'][iwan]['m']*2)==1):
                                control['impurity_wan'][ip][10]=wan_hmat['basis'][iwan]['ind']
                            elif (int(wan_hmat['basis'][iwan]['m']*2)==3):
                                control['impurity_wan'][ip][11]=wan_hmat['basis'][iwan]['ind']
                            elif (int(wan_hmat['basis'][iwan]['m']*2)==5):
                                control['impurity_wan'][ip][12]=wan_hmat['basis'][iwan]['ind']
                            elif (int(wan_hmat['basis'][iwan]['m']*2)==7):
                                control['impurity_wan'][ip][13]=wan_hmat['basis'][iwan]['ind']
                if (control['impurity_wan'][ip].count(0) !=0):
                    print('something wrong in find_impurity_wan', file=control['h_log'],flush=True)
                    sys.exit()
        else:
            if   (control['impurity_problem'][ip][1].lower()=='s'):
                control['impurity_wan'].append([0]*1)
                for iwan in range(num_wann):
                    if ((wan_hmat['basis'][iwan]['atom']==control['impurity_problem'][ip][0]) and (wan_hmat['basis'][iwan]['l']==0)):
                        if (wan_hmat['basis'][iwan]['m']==-0):
                            control['impurity_wan'][ip][0]=wan_hmat['basis'][iwan]['ind']
                if (control['impurity_wan'][ip].count(0) !=0):
                    print('something wrong in find_impurity_wan', file=control['h_log'],flush=True)
                    sys.exit()    
            elif (control['impurity_problem'][ip][1].lower()=='p'):
                control['impurity_wan'].append([0]*3)
                for iwan in range(num_wann):
                    if ((wan_hmat['basis'][iwan]['atom']==control['impurity_problem'][ip][0]) and (wan_hmat['basis'][iwan]['l']==1)):
                        if (wan_hmat['basis'][iwan]['m']==-1):
                            control['impurity_wan'][ip][0]=wan_hmat['basis'][iwan]['ind']
                        elif (wan_hmat['basis'][iwan]['m']==-0):
                            control['impurity_wan'][ip][1]=wan_hmat['basis'][iwan]['ind']
                        elif (wan_hmat['basis'][iwan]['m']==1):
                            control['impurity_wan'][ip][2]=wan_hmat['basis'][iwan]['ind']
                if (control['impurity_wan'][ip].count(0) !=0):
                    print('something wrong in find_impurity_wan', file=control['h_log'],flush=True)
                    sys.exit()    
            elif (control['impurity_problem'][ip][1].lower()=='d'):
                control['impurity_wan'].append([0]*5)
                for iwan in range(num_wann):
                    if ((wan_hmat['basis'][iwan]['atom']==control['impurity_problem'][ip][0]) and (wan_hmat['basis'][iwan]['l']==2)):
                        if (wan_hmat['basis'][iwan]['m']==-2):
                            control['impurity_wan'][ip][0]=wan_hmat['basis'][iwan]['ind']
                        elif (wan_hmat['basis'][iwan]['m']==-1):
                            control['impurity_wan'][ip][1]=wan_hmat['basis'][iwan]['ind']
                        elif (wan_hmat['basis'][iwan]['m']==-0):
                            control['impurity_wan'][ip][2]=wan_hmat['basis'][iwan]['ind']
                        elif (wan_hmat['basis'][iwan]['m']==1):
                            control['impurity_wan'][ip][3]=wan_hmat['basis'][iwan]['ind']
                        elif (wan_hmat['basis'][iwan]['m']==2):
                            control['impurity_wan'][ip][4]=wan_hmat['basis'][iwan]['ind']
                if (control['impurity_wan'][ip].count(0) !=0):
                    print('something wrong in find_impurity_wan', file=control['h_log'],flush=True)
                    sys.exit()    
            elif (control['impurity_problem'][ip][1].lower()=='f'):
                control['impurity_wan'].append([0]*7)
                for iwan in range(num_wann):
                    if ((wan_hmat['basis'][iwan]['atom']==control['impurity_problem'][ip][0]) and (wan_hmat['basis'][iwan]['l']==3)):
                        if (wan_hmat['basis'][iwan]['m']==-3):
                            control['impurity_wan'][ip][0]=wan_hmat['basis'][iwan]['ind']
                        elif (wan_hmat['basis'][iwan]['m']==-2):
                            control['impurity_wan'][ip][1]=wan_hmat['basis'][iwan]['ind']
                        elif (wan_hmat['basis'][iwan]['m']==-1):
                            control['impurity_wan'][ip][2]=wan_hmat['basis'][iwan]['ind']
                        elif (wan_hmat['basis'][iwan]['m']==-0):
                            control['impurity_wan'][ip][3]=wan_hmat['basis'][iwan]['ind']
                        elif (wan_hmat['basis'][iwan]['m']==1):
                            control['impurity_wan'][ip][4]=wan_hmat['basis'][iwan]['ind']
                        elif (wan_hmat['basis'][iwan]['m']==2):
                            control['impurity_wan'][ip][5]=wan_hmat['basis'][iwan]['ind']
                        elif (wan_hmat['basis'][iwan]['m']==3):
                            control['impurity_wan'][ip][6]=wan_hmat['basis'][iwan]['ind']
                if (control['impurity_wan'][ip].count(0) !=0):
                    print('something wrong in find_impurity_wan', file=control['h_log'],flush=True)
                    sys.exit()    

    return None




def initial_file_directory_setup(control):

    directory_setup(control)
    if (control['method'] == 'lda+dmft'):
        print('iter_num_impurity', control['iter_num_impurity'], '  max_iter_num_impurity', control['max_iter_num_impurity'], file=control['h_log'],flush=True)
        print('iter_num_outer', control['iter_num_outer'], '  max_iter_num_outer', control['max_iter_num_outer'],             file=control['h_log'],flush=True)
    elif (control['method'] == 'lqsgw+dmft'):
        print('iter_num_impurity', control['iter_num_impurity'],         file=control['h_log'],flush=True)
        print('max_iter_num_impurity', control['max_iter_num_impurity'], file=control['h_log'],flush=True)

    return None




def find_place_to_restart(control):

    if (control['method']=='lqsgw+dmft'):
        control['conv_table']=read_convergence_table(control)
        # print(control['conv_table'], file=control['h_log'],flush=True)

        if (len(control['conv_table'])>0):

            n_imp_problem=amax(control['impurity_problem_equivalence'])

            last_step=control['conv_table'][-1][0].strip().split('_')[0]
            last_imp_iter=control['conv_table'][-1][1].strip()
            if (len(control['conv_table'][-1][0].strip().split('_')) > 1):
                last_imp=control['conv_table'][-1][0].strip().split('_')[1]
                print(last_step, last_imp, last_imp_iter, file=control['h_log'],flush=True)
            else:
                print(last_step, last_imp_iter, file=control['h_log'],flush=True)
            if last_step == 'wannier':
                control['do_wannier']=False
                control['do_coulomb']=True
                control['do_dc']=True
                control['iter_num_impurity']=1
            elif last_step == 'coulomb':
                control['do_wannier']=False
                control['do_coulomb']=False
                control['do_dc']=True
                control['iter_num_impurity']=1
            elif last_step == 'dc':
                if (int(last_imp) == n_imp_problem):
                    control['do_wannier']=False
                    control['do_coulomb']=False
                    control['do_dc']=False
                    control['iter_num_impurity']=1
                else:
                    control['do_wannier']=False
                    control['do_coulomb']=False
                    control['do_dc']=True
                    control['iter_num_impurity']=1
                    for ii in range(int(last_imp)):
                        control['conv_table'].pop(-1)
            elif (last_step == 'delta'):
                control['do_wannier']=False
                control['do_coulomb']=False
                control['do_dc']=False
                control['iter_num_impurity']=int(last_imp_iter)
                control['conv_table'].pop(-1)
            elif (last_step == 'impurity'):
                if (int(last_imp) == n_imp_problem):
                    control['do_wannier']=False
                    control['do_coulomb']=False
                    control['do_dc']=False
                    control['iter_num_impurity']=int(last_imp_iter)+1
                else:
                    control['do_wannier']=False
                    control['do_coulomb']=False
                    control['do_dc']=True
                    control['iter_num_impurity']=int(last_imp_iter)
                    for ii in range(int(last_imp)):
                        control['conv_table'].pop(-1)
            else:
                control['do_wannier']=True
                control['do_coulomb']=True
                control['do_dc']=True
                control['iter_num_impurity']=1
        else:
            control['do_wannier']=True
            control['do_coulomb']=True
            control['do_dc']=True
            control['iter_num_impurity']=1
    elif (control['method']=='lda+dmft'):
        control['conv_table']=read_convergence_table(control)
        if (len(control['conv_table'])>0):
            for ii in range(np.shape(control['conv_table'])[0]):
                if control['conv_table'][ii][0].strip()=='dft':
                    linecnt=ii
                    control['iter_num_outer']=int(control['conv_table'][ii][1])
            for ii in range(linecnt, np.shape(control['conv_table'])[0]):
                control['conv_table'].pop(-1)
    return None

# def find_iter_num_for_restart(control):

#     if (control['restart']):
#         line_count=sum(1 for line in open(control['top_dir']+'/convergence.log'))        
#         if (line_count <=1):
#             if (control['method']=='lda+dmft'):
#                 iter_num_outer=1
#             elif (control['method']=='lqsgw+dmft'):
#                 iter_num_impurity=1                                
#         else:
#             if (control['method']=='lda+dmft'):
#                 iter_num_outer=1                
#                 ff=open(control['top_dir']+'/convergence.log', 'r')
#                 firstline=ff.readline()                
#                 for line in ff:
#                     temp=line.split()
#                     if (temp[0] == 'dft'):
#                         iter_num_outer=int(temp[1])
#                 ff.close()
#             elif (control['method']=='lqsgw+dmft'):
#                 iter_num_impurity=1                                               
#                 ff=open(control['top_dir']+'/convergence.log', 'r')
#                 firstline=ff.readline()                            
#                 for line in ff:
#                     temp=line.split()
#                     temp1=temp[0]
#                     if (temp1 == 'impurity'):
#                         iter_num_impurity=int(temp[2])
#                 ff.close()
#     else:
#         if (control['method']=='lda+dmft'):
#             iter_num_outer=1
#         elif (control['method']=='lqsgw+dmft'):
#             iter_num_impurity=1

#     if (control['method']=='lda+dmft'):
#         return iter_num_outer
#     elif (control['method']=='lqsgw+dmft'):
#         return iter_num_impurity





def initial_lattice_directory_setup(control):
    os.chdir(control['lattice_directory'])
    files = glob.iglob(control['initial_lattice_dir']+"/*.rst")
    for filename in files:
        shutil.copy(filename, './')
    files = glob.iglob(control['initial_lattice_dir']+"/*el_density")
    for filename in files:
        shutil.copy(filename, './')    
    if os.path.exists(control['initial_lattice_dir']+'/kpath'):
        shutil.copy(control['initial_lattice_dir']+'/kpath', './')
    if os.path.exists(control['initial_lattice_dir']+'/ini'):
        shutil.copy(control['initial_lattice_dir']+'/ini', './')

    iter_string='_'+str(control['iter_num_outer'])

    shutil.copy(control['initial_lattice_dir']+'/'+control['allfile']+'.out', control['allfile']+iter_string+'.out')

    print("initial dft directory setup done", file=control['h_log'],flush=True)
    os.chdir(control['top_dir'])
    return None

def create_comwann_ini(control, wan_hmat):
    f=open('comwann.ini','w')

    if (control['method']=='lda+dmft'):
        f.write(control['lattice_directory']+'\n')
        f.write('dft\n')
    elif (control['method']=='lqsgw+dmft'):
        f.write(control['initial_lattice_dir']+'\n')
        f.write('qp\n')
    f.write(str(wan_hmat['dis_win_max'])+'\n')
    f.write(str(wan_hmat['dis_win_min'])+'\n')
    f.write(str(wan_hmat['froz_win_max'])+'\n')
    f.write(str(wan_hmat['froz_win_min'])+'\n')

    f.write(str(wan_hmat['num_iter'])+'\n')
    f.write(str(wan_hmat['dis_num_iter'])+'\n')
    f.write('0\n')
    f.write(str(wan_hmat['cut_low'])+'\n')
    f.write(str(wan_hmat['cut_froz'])+'\n')
    f.write(str(wan_hmat['cut_total'])+'\n')    
    f.write(str(wan_hmat['rmode'])+'\n')
    f.write(str(wan_hmat['radfac'])+'\n')    

    f.close()

def create_comcoulomb_ini(control):
    f=open('comcoulomb.ini','w')

    f.write(control['initial_lattice_dir']+'\n')
    f.write(control['wannier_directory']+'\n')

    f.write(str(control['nproc_tau_coulomb'])+'\n')
    f.write(str(control['nproc_k_coulomb'])+'\n')    
    f.write(str(control['proj_win_min'])+'\n')
    f.write(str(control['proj_win_max'])+'\n')
    f.write('F\n')
    f.write(control['u_mode']+'\n')
    nimp_orb=0
    natom=len(control['impurity_wan'])    
    for ii in range(natom):
        nimp_orb=nimp_orb+len(control['impurity_wan'][ii])
    f.write(str(nimp_orb)+'\n')
    for iatom in range(natom):
        f.write(' '.join(map(str,control['impurity_wan'][iatom]))+' ')
    f.write('\n')
    f.write('1\n')
    f.write('F\n')
    f.write('3.0\n')
    f.write('F\n')            

    f.close()    


# def create_wannier_inip(wan_hmat):
#     # in the wannier directory    
#     g=open('wannier.inip', 'w')
#     num_wann=np.shape(wan_hmat['basis'])[0]
#     g.write(str(num_wann)+'\n')
#     for ii in range(num_wann):
#         if (control['spin_orbit']==False):            
#             tempstr=[wan_hmat['basis'][ii]['atom'], wan_hmat['basis'][ii]['l'], wan_hmat['basis'][ii]['m'], wan_hmat['basis'][ii]['xaxis'][0], wan_hmat['basis'][ii]['xaxis'][1], wan_hmat['basis'][ii]['xaxis'][2], wan_hmat['basis'][ii]['zaxis'][0], wan_hmat['basis'][ii]['zaxis'][1], wan_hmat['basis'][ii]['zaxis'][2]]
#         else:
#             tempstr=[wan_hmat['basis'][ii]['atom'], wan_hmat['basis'][ii]['l'], wan_hmat['basis'][ii]['i'], wan_hmat['basis'][ii]['m'], wan_hmat['basis'][ii]['xaxis'][0], wan_hmat['basis'][ii]['xaxis'][1], wan_hmat['basis'][ii]['xaxis'][2], wan_hmat['basis'][ii]['zaxis'][0], wan_hmat['basis'][ii]['zaxis'][1], wan_hmat['basis'][ii]['zaxis'][2]]
#         g.write(' '.join(map(str, tempstr))+'\n')            
#     g.close()

#     return None


def read_wan_hmat_basis(control):
    # in the wannier directory

    inip=loadtxt(control['wannier_directory']+'/wannier.inip')
    basis_info=[]

    if (control['spin_orbit']):
        for ii in range(np.shape(inip)[0]):
            basis_info.append({'atom':int(inip[ii,0]), 'l':int(inip[ii,1]), 'i':inip[ii,2],'m':inip[ii,3],'xaxis':inip[ii,4:7],'zaxis':inip[ii,7:10], 'ind':ii+1})
    else:
        for ii in range(np.shape(inip)[0]):
            basis_info.append({'atom':int(inip[ii,0]), 'l':int(inip[ii,1]), 'm':int(inip[ii,2]),'xaxis':inip[ii,3:6],'zaxis':inip[ii,6:9], 'ind':ii+1})

    print(basis_info,                                      file=control['h_log'],flush=True)
    print('reading wannier.inip to get basis information', file=control['h_log'],flush=True)
    return basis_info

def check_key_in_string(key,dictionary,control):
    if (key not in dictionary):
        print('missing \''+key+'\' in '+dictionary['name'], file=control['h_log'],flush=True)
        sys.exit()
    return None


def overwrite_key_in_string(key,dictionary,dictionaryname,value,h_log):
    if (key in dictionary):
        print('\''+key+'\' in '+dictionaryname+' is overwritten', file=control['h_log'],flush=True)
    return value

# def dft_rst_file_check():
#     check_for_files('*acc_core_dft.rst', h_log)
#     check_for_files('*chemical_potential_dft.rst', h_log)
#     check_for_files('*cor_norm_dft.rst', h_log)
#     check_for_files('*dfi_dft.rst', h_log)
#     check_for_files('*dfidot2_dft.rst', h_log)
#     check_for_files('*dfidot_dft.rst', h_log)
#     check_for_files('*e_bnd_dft.rst', h_log)
#     check_for_files('*e_core_dft.rst', h_log)
#     check_for_files('*el_density_dft.rst', h_log)
#     check_for_files('*eny_dft.rst', h_log)
#     check_for_files('*etot_dft.rst', h_log)
#     check_for_files('*ev_bnd_*_dft.rst', h_log)
#     check_for_files('*ffsmt_dft.rst', h_log)
#     check_for_files('*fi_dft.rst', h_log)
#     check_for_files('*fidot2_dft.rst', h_log)
#     check_for_files('*fidot_dft.rst', h_log)
#     check_for_files('*g_full_00_*_dft.rst', h_log)
#     check_for_files('*g_loc_0_dft.rst', h_log)
#     check_for_files('*gfun_dft.rst', h_log)
#     check_for_files('*gfun_old_dft.rst', h_log)
#     check_for_files('*gfund_dft.rst', h_log)
#     check_for_files('*gfund_old_dft.rst', h_log)
#     check_for_files('*n_bnd_dft.rst', h_log)
#     check_for_files('*p_f_dft.rst', h_log)
#     check_for_files('*pcor_dft.rst', h_log)
#     check_for_files('*pcor_old_dft.rst', h_log)
#     check_for_files('*pd2_f_dft.rst', h_log)
#     check_for_files('*pd_f_dft.rst', h_log)
#     check_for_files('*ptnl_dft.rst', h_log)
#     check_for_files('*q_f_dft.rst', h_log)
#     check_for_files('*qcor_dft.rst', h_log)
#     check_for_files('*qcor_old_dft.rst', h_log)
#     check_for_files('*qd2_f_dft.rst', h_log)
#     check_for_files('*qd_f_dft.rst', h_log)
#     check_for_files('*restart_ubi.rst', h_log)
#     check_for_files('*ro_core_dft.rst', h_log)
#     check_for_files('*v_intr_h_dft.rst', h_log)
#     check_for_files('*v_intr_xc_dft.rst', h_log)
#     check_for_files('*v_mt_h_dft.rst', h_log)
#     check_for_files('*v_mt_xc_dft.rst', h_log)
#     check_for_files('*z_bnd_*_dft.rst', h_log)
#     return None


# def string_addwhitespace(string, stringsize):
#     stringout=string
#     if stringsize > len(string):
#         stringout=string+' '*(stringsize-len(string))
#     return stringout


def find_all_in_string(str, ch):
    for i, ltr in enumerate(str):
        if ltr == ch:
            yield i


def read_convergence_table(control):

    if os.path.exists(control['top_dir']+'/convergence.log'):
        with open(control['top_dir']+'/convergence.log', 'r') as logfile:
            tmp=logfile.readlines()
        nstep=len(tmp)-2
        if (nstep>0):
            endind=list(find_all_in_string(tmp[1],' '))[::2]+[len(tmp[1])-1]
            startind=[0]+(array(list(find_all_in_string(tmp[1],' '))[1::2])+1).tolist()
            ncolumn=len(endind)
            f=open('./convergence.log', 'r')
            f.readline()
            f.readline()
            convergence_table=[]
            for lines in f:
                eachline=[]
                for ii in range(ncolumn):
                    eachline.append(lines.rstrip()[startind[ii]:endind[ii]])
                if (len(eachline[0])>0):
                    convergence_table.append(eachline)
            f.close()
        else:
            convergence_table=[]
    else:
        convergence_table=[]    
    return convergence_table


def generate_initial_self_energy(control,imp):
    os.chdir(control['impurity_directory'])    
    if ('initial_self_energy' in control):
        shutil.copy(control['initial_self_energy'], './sig.dat')
        if ('initial_impurity_dir' in control):
            initial_impurity_dirname=os.path.abspath(os.path.dirname(control['initial_impurity_dir']))
            directories = glob.glob(initial_impurity_dirname+"/*/")
            for directory_name in directories:
                dest_dir=directory_name.split('/')[-2]
                files = glob.iglob(os.path.abspath(directory_name)+"/config*")
                for filename in files:
                    shutil.copy(filename, control['impurity_directory']+'/'+dest_dir)
    else:
        dc=np.loadtxt(control['dc_directory']+'/dc.dat')
        beta=imp['beta']
        n_omega=control['n_omega']
        omega=control['omega']


        sig_table=[]
        for jj in range(control['n_omega']):
            sig_omega=[control['omega'][jj]]+dc.tolist()
            sig_table.append(sig_omega)
        with open('./sig.dat', 'w') as outputfile:
            outputfile.write(tabulate(sig_table, headers=control['sig_header'], floatfmt=".12f", numalign="right",  tablefmt="plain"))
    if (control['method']=='lqsgw+dmft'):
        iter_string='_0'
    elif (control['method']=='lda+dmft'):
        iter_string='_'+str(control['iter_num_outer'])+'_0'

    labeling_file('./sig.dat', iter_string)
    print('initial_self_energy generation done', file=control['h_log'],flush=True)
    os.chdir(control['top_dir'])
    return None

def prepare_initial_ef(control):
    os.chdir(control['lowh_directory'])
    f=open('ef.dat','w')
    f.write('0.0\n')
    f.close()
    os.chdir(control['top_dir'])
    return None    


def delta_postprocessing(control,imp):

    write_transformation_matrix(control,control['lowh_directory']+'/local_spectral_matrix_ef.dat')
    cal_projected_mean_field_diagonal(control,imp)
    cal_dc_diagonal(control)
    cal_e_imp_diagonal(control)
    delta_causality=cal_hyb_diagonal(control,imp)

    if (delta_causality ==0):
        print('delta causality broken', file=control['h_log'],flush=True)
        sys.exit()    

    return delta_causality


def cal_dc_diagonal(control):

    os.chdir(control['dc_directory'])
    dc_mat=read_impurity_mat_static(control,control['dc_directory']+'/dc_mat.dat')

    h=open('./dc.dat', 'w')

    for ii in sorted(set(control['impurity_problem_equivalence'])):
        dc_vec=imp_from_mat_to_array(dc_mat[str(ii)],imp[str(abs(ii))]['impurity_matrix'])
        for jj in range(len(dc_vec)):
            h.write(str(real(dc_vec[jj]))+'   '+str(imag(dc_vec[jj]))+'    ')
    h.close()

    if (control['method']=='lqsgw+dmft'):
        iter_string='_'+str(control['iter_num_impurity'])
    elif (control['method']=='lda+dmft'):
        iter_string='_'+str(control['iter_num_outer'])+'_'+str(control['iter_num_impurity'])

    labeling_file('./dc.dat', iter_string)

    print('dc.dat generation done', file=control['h_log'],flush=True)
    os.chdir(control['top_dir'])

    return None


def cal_dc_diagonal_new(control):

    os.chdir(control['dc_directory'])
    dc_mat=read_impurity_mat_static(control,control['dc_directory']+'/dc_mat.dat')

    h=open('./dc.dat', 'w')

    for ii in sorted(set(control['impurity_problem_equivalence'])):
        dc_vec=imp_from_mat_to_array(dc_mat[str(ii)],imp[str(abs(ii))]['impurity_matrix'])
        for jj in range(len(dc_vec)):
            h.write(str(real(dc_vec[jj]))+'   '+str(imag(dc_vec[jj]))+'    ')
    h.close()

    if (control['method']=='lqsgw+dmft'):
        iter_string='_'+str(control['iter_num_impurity'])
    elif (control['method']=='lda+dmft'):
        iter_string='_'+str(control['iter_num_outer'])+'_'+str(control['iter_num_impurity'])

    labeling_file('./dc.dat', iter_string)

    print('dc.dat generation done', file=control['h_log'],flush=True)
    os.chdir(control['top_dir'])

    return None

def cal_zinv_m1_diagonal(control):

    os.chdir(control['dc_directory'])
    zinv_m1_mat=read_impurity_mat_static(control,control['dc_directory']+'/zinv_m1_mat.dat')

    h=open('./zinv_m1.dat', 'w')

    for ii in sorted(set(control['impurity_problem_equivalence'])):
        zinv_m1_vec=imp_from_mat_to_array(zinv_m1_mat[str(ii)],imp[str(abs(ii))]['impurity_matrix'])        
        for jj in range(len(zinv_m1_vec)):
            h.write(str(real(zinv_m1_vec[jj]))+'   '+str(imag(zinv_m1_vec[jj]))+'    ')
    h.close()

    if (control['method']=='lqsgw+dmft'):
        iter_string='_'+str(control['iter_num_impurity'])
    elif (control['method']=='lda+dmft'):
        iter_string='_'+str(control['iter_num_outer'])+'_'+str(control['iter_num_impurity'])

    labeling_file('./zinv_m1.dat', iter_string)

    print('zinv_m1.dat generation done', file=control['h_log'],flush=True)
    os.chdir(control['top_dir'])

    return None

def vec_from_mat_dynamic(mat,trans):
    vec=zeros(np.shape(mat, 0), np.shape(mat, 1))
    for ii in range(np.shape(mat, 0)):
        vec[ii,:]=np.diag(dot(np.transpose(np.conj(trans)), np.dot(mat[ii,:,:], trans)))
    return vec


def prepare_impurity_solver(control,wan_hmat,imp):

    # cal_trans_from_patrick(control, imp)
    delta=array_impurity_dynamic(control,imp,control['lowh_directory']+'/delta.dat')

    write_json_all(control,imp,delta,'hyb.json')

    e_imp=generate_mat_from_array_impurity_static(control,imp,control['lowh_directory']+'/e_imp.dat')
    trans_basis=read_impurity_mat_static(control,control['lowh_directory']+'/trans_basis.dat')

    for key, value in imp.items(): 
        if (not (isinstance(imp[key], dict))):
            continue
        nimp_orb=len(imp[key]['impurity_matrix'])
        os.chdir(control['impurity_directory']+'/'+key)
        if (control['spin_orbit']):
            ndim=nimp_orb
            e_imp_key=np.zeros((ndim, ndim))
            trans_key=np.zeros((ndim, ndim))
            equivalence_key=np.zeros((ndim,ndim),dtype=integer)
            e_imp_key=real(e_imp[key])
            trans_key=real(trans_basis[key])
            # equivalence_key=array([[(lambda ii: str(ii) if str(ii)!='0' else '')(ii) for ii in row] for row in imp[key]['impurity_matrix']])
            equivalence_key=list(map(lambda row: list(map(lambda x: str(x) if x!='0' else '', list(map(str, row)))), imp[key]['impurity_matrix']))
        else:
            ndim=nimp_orb*2
            e_imp_key=np.zeros((ndim, ndim))
            trans_key=np.zeros((ndim, ndim))
            equivalence_key=np.zeros((ndim,ndim),dtype='<U2')
            equivalence_key_int_mat=array(imp[key]['impurity_matrix'])
            print(equivalence_key_int_mat, file=control['h_log'],flush=True)
            if (imp[key]['para']):
                mkey=key
                shiftval=0
            else:
                mkey=str(-int(key))
                shiftval=amax(equivalence_key_int_mat)
            print(mkey, shiftval, file=control['h_log'],flush=True)
            #
            # On the next line ii>0 evaluates to 1 if ii>0 and evaluates to 0 otherwise
            # equivalence_mkey_int_mat=equivalence_key_int_mat+shiftval*array([[(lambda ii: ii>0)(ii) for ii in row] for row in equivalence_key_int_mat])
            # equivalence_mkey_int_mat=equivalence_key_int_mat+shiftval*array(map(lambda row: map(int,row), equivalence_key_int_mat>0))
            equivalence_mkey_int_mat=equivalence_key_int_mat+shiftval*(equivalence_key_int_mat>0)


            e_imp_key[0:nimp_orb,0:nimp_orb]=real(e_imp[key])
            e_imp_key[nimp_orb:(2*nimp_orb),nimp_orb:(2*nimp_orb)]=real(e_imp[mkey])
            trans_key[0:nimp_orb,0:nimp_orb]=real(trans_basis[key])
            trans_key[nimp_orb:(2*nimp_orb),nimp_orb:(2*nimp_orb)]=real(trans_basis[mkey])
            # equivalence_key[0:nimp_orb,0:nimp_orb]=[[(lambda ii: str(ii) if str(ii)!='0' else '')(ii) for ii in row] for row in equivalence_key_int_mat]    
            equivalence_key[0:nimp_orb,0:nimp_orb]=list(map(lambda row: list(map(lambda x: str(x) if x!='0' else '', list(map(str, row)))), equivalence_key_int_mat))
            # equivalence_key[nimp_orb:(2*nimp_orb),nimp_orb:(2*nimp_orb)]=[[(lambda ii: str(ii) if str(ii)!='0' else '')(ii) for ii in row] for row in equivalence_key_int_mat]
            equivalence_key[nimp_orb:(2*nimp_orb),nimp_orb:(2*nimp_orb)]=list(map(lambda row: list(map(lambda x: str(x) if x!='0' else '', list(map(str, row)))), equivalence_mkey_int_mat))

        write_params_json(control,imp[key],e_imp_key,trans_key,equivalence_key,imp['beta'])
        if (control['method']=='lqsgw+dmft'):
            write_dynamical_f0_json(imp[key])

    os.chdir(control['top_dir'])
    return None


def run_impurity_solver(control,imp):
    green={}
    sigma_bare={}
    sigma={}
    sigma_to_delta={}
    for key, value in imp.items():
        if (not (isinstance(imp[key], dict))):
            continue
        os.chdir(control['impurity_directory']+'/'+key)
        solve_impurity_patrick(control)
        measure_impurity_patrick(control)
        green[key], sigma_bare[key], sigma[key], sigma_to_delta[key]=impurity_postprocessing(control, imp, key)

    os.chdir(control['impurity_directory'])


    green_table=[]
    sigma_table=[]
    sigma_to_delta_table=[]
    sigma_bare_table=[]        
    for jj in range(control['n_omega']):
        green_omega=[control['omega'][jj]]
        sigma_omega=[control['omega'][jj]]
        sigma_to_delta_omega=[control['omega'][jj]]
        sigma_bare_omega=[control['omega'][jj]]
        for ii in sorted(set(control['impurity_problem_equivalence'])):
            n_iio=amax(imp[str(abs(ii))]['impurity_matrix'])
            for kk in range(n_iio):
                if (ii<0):
                    pp=kk+n_iio
                else:
                    pp=kk
                green_omega=green_omega+[real(green[str(abs(ii))][jj,pp]),imag(green[str(abs(ii))][jj,pp])]
                sigma_omega=sigma_omega+[real(sigma[str(abs(ii))][jj,pp]),imag(sigma[str(abs(ii))][jj,pp])]
                sigma_to_delta_omega=sigma_to_delta_omega+[real(sigma_to_delta[str(abs(ii))][jj,pp]),imag(sigma_to_delta[str(abs(ii))][jj,pp])]
                sigma_bare_omega=sigma_bare_omega+[real(sigma_bare[str(abs(ii))][jj,pp]),imag(sigma_bare[str(abs(ii))][jj,pp])]
        green_table.append(green_omega)
        sigma_table.append(sigma_omega)
        sigma_to_delta_table.append(sigma_to_delta_omega)
        sigma_bare_table.append(sigma_bare_omega)

    with open('./gimp.dat', 'w') as outputfile:
        outputfile.write(tabulate(green_table, headers=control['sig_header'], floatfmt=".12f", numalign="right",  tablefmt="plain"))
    with open('./sig_bare.dat', 'w') as outputfile:
        outputfile.write(tabulate(sigma_bare_table, headers=control['sig_header'], floatfmt=".12f", numalign="right",  tablefmt="plain"))
    with open('./sig_smth.dat', 'w') as outputfile:
        outputfile.write(tabulate(sigma_table, headers=control['sig_header'], floatfmt=".12f", numalign="right",  tablefmt="plain"))
    with open('./sig.dat', 'w') as outputfile:
        outputfile.write(tabulate(sigma_to_delta_table, headers=control['sig_header'], floatfmt=".12f", numalign="right",  tablefmt="plain"))

    shutil.copy('./sig.dat', control['top_dir'])

    if (control['method']=='lqsgw+dmft'):
        iter_string='_'+str(control['iter_num_impurity'])
    elif (control['method']=='lda+dmft'):
        iter_string='_'+str(control['iter_num_outer'])+'_'+str(control['iter_num_impurity'])    

    labeling_file('./gimp.dat',iter_string)
    labeling_file('./sig_bare.dat',iter_string)
    labeling_file('./sig_smth.dat',iter_string)
    labeling_file('./sig.dat',iter_string)

    os.chdir(control['top_dir'])


def generate_mat_from_array_impurity_dynamic(control,imp, filename):

    os.chdir(control['impurity_directory'])

    dat=loadtxt(filename)

    start_array={}
    end_array={}

    last_index=1

    for ii in sorted(set(control['impurity_problem_equivalence'])):
        n_iio=amax(imp[str(abs(ii))]['impurity_matrix'])
        start_array[ii]=last_index
        end_array[ii]=last_index+2*n_iio
        last_index=last_index+2*n_iio
    print(start_array)
    print(end_array)

    matout={}

    for ii in sorted(set(control['impurity_problem_equivalence'])):
        nimp_orb=len(imp[str(abs(ii))]['impurity_matrix'])
        tempmat=np.zeros((control['n_omega'],nimp_orb,nimp_orb), dtype='complex')
        for iomega in range(control['n_omega']):
            tempmat2=dat[iomega,start_array[ii]:end_array[ii]]
            tempmat[iomega,:,:]=imp_from_array_to_mat(tempmat2[0::2]+tempmat2[1::2]*1j,imp[str(abs(ii))]['impurity_matrix'])
        matout[str(ii)]=tempmat
    return matout


def generate_mat_from_array_impurity_static(control,imp, filename):

    os.chdir(control['impurity_directory'])

    dat=loadtxt(filename)

    start_array={}
    end_array={}

    last_index=0

    for ii in sorted(set(control['impurity_problem_equivalence'])):
        n_iio=amax(imp[str(abs(ii))]['impurity_matrix'])
        start_array[ii]=last_index
        end_array[ii]=last_index+2*n_iio
        last_index=last_index+2*n_iio
    print(start_array)
    print(end_array)

    matout={}
    for ii in sorted(set(control['impurity_problem_equivalence'])):
        tempmat2=dat[start_array[ii]:end_array[ii]]
        matout[str(ii)]=imp_from_array_to_mat(tempmat2[0::2]+tempmat2[1::2]*1j,imp[str(abs(ii))]['impurity_matrix'])
    return matout


def array_impurity_static(control,imp, filename):

    os.chdir(control['impurity_directory'])

    dat=loadtxt(filename)

    start_array={}
    end_array={}

    last_index=0

    for ii in sorted(set(control['impurity_problem_equivalence'])):
        n_iio=amax(imp[str(abs(ii))]['impurity_matrix'])
        start_array[ii]=last_index
        end_array[ii]=last_index+2*n_iio
        last_index=last_index+2*n_iio
    print(start_array)
    print(end_array)

    matout={}
    for ii in sorted(set(control['impurity_problem_equivalence'])):
        tempmat2=dat[start_array[ii]:end_array[ii]]
        matout[str(ii)]=tempmat2[0::2]+tempmat2[1::2]*1j
    return matout


def array_impurity_dynamic(control,imp, filename):

    os.chdir(control['impurity_directory'])

    dat=loadtxt(filename)

    start_array={}
    end_array={}

    last_index=1

    for ii in sorted(set(control['impurity_problem_equivalence'])):
        n_iio=amax(imp[str(abs(ii))]['impurity_matrix'])
        start_array[ii]=last_index
        end_array[ii]=last_index+2*n_iio
        last_index=last_index+2*n_iio
    print(start_array)
    print(end_array)

    matout={}

    for ii in sorted(set(control['impurity_problem_equivalence'])):
        n_iio=amax(imp[str(abs(ii))]['impurity_matrix'])
        tempmat=np.zeros((control['n_omega'],n_iio), dtype='complex')
        for iomega in range(control['n_omega']):
            tempmat2=dat[iomega,start_array[ii]:end_array[ii]]
            tempmat[iomega,:]=tempmat2[0::2]+tempmat2[1::2]*1j
        matout[str(ii)]=tempmat
    return matout

def cal_projected_mean_field_diagonal(control,imp):

    os.chdir(control['lowh_directory'])
    hmat=read_impurity_mat_static(control,control['lowh_directory']+'/e_projected_mat.dat')

    h=open('./projected_eig.dat', 'w')


    for ii in sorted(set(control['impurity_problem_equivalence'])):    
        h_vec=imp_from_mat_to_array(hmat[str(ii)],imp[str(abs(ii))]['impurity_matrix'])

        for jj in range(len(h_vec)):
            h.write(str(real(h_vec[jj]))+'   '+str(imag(h_vec[jj]))+'    ')
    h.close()

    if (control['method']=='lqsgw+dmft'):
        iter_string='_'+str(control['iter_num_impurity'])
    elif (control['method']=='lda+dmft'):
        iter_string='_'+str(control['iter_num_outer'])+'_'+str(control['iter_num_impurity'])

    labeling_file('./projected_eig.dat', iter_string)

    print('projected_eig.dat generation done', file=control['h_log'],flush=True)
    os.chdir(control['top_dir'])
    return None

def cal_e_imp_diagonal(control):
    os.chdir(control['lowh_directory'])
    eig=np.loadtxt('projected_eig.dat')
    dc=np.loadtxt(control['dc_directory']+'/dc.dat')

    f=open('e_imp.dat', 'w')
    f.write("  ".join(map(str, eig-dc))+'\n')
    f.close()

    if (control['method']=='lqsgw+dmft'):
        iter_string='_'+str(control['iter_num_impurity'])
    elif (control['method']=='lda+dmft'):
        iter_string='_'+str(control['iter_num_outer'])+'_'+str(control['iter_num_impurity'])

    labeling_file('./e_imp.dat', iter_string)

    print('e_imp.dat generation done', file=control['h_log'],flush=True)
    os.chdir(control['top_dir'])
    return None


def imp_from_array_to_mat(vecin,equivalence_mat):
    nimp_orb=len(equivalence_mat)
    matout=np.zeros((nimp_orb, nimp_orb), dtype='complex')
    for ii in range(nimp_orb):    
        for jj in range(nimp_orb):
            if (equivalence_mat[ii,jj]!=0):
                matout[ii,jj]=vecin[equivalence_mat[ii,jj]-1]
    return matout



def imp_from_mat_to_array(matin,equivalence_mat):
    n_iio=amax(equivalence_mat)
    vecout=np.zeros(n_iio, dtype='complex')
    degen_vec=np.zeros(n_iio, dtype='int')
    nimp_orb=len(matin)
    # print(nimp_orb)
    # print(equivalence_mat)
    # print(type(equivalence_mat))
    # print(matin)
    # print(type(matin))

    for ii in range(nimp_orb):    
        for jj in range(nimp_orb):
            print(ii, jj)
            if (equivalence_mat[ii,jj]!=0):
                ind=equivalence_mat[jj,jj]-1
                vecout[ind]=vecout[ind]+matin[ii,jj]
                degen_vec[ind]=degen_vec[ind]+1
    vecout=vecout/(degen_vec*1.0)

    return vecout

# def read_trans_basis(control,filename):
#     trans_basis={}
#     g=open(filename, 'r')

#     for ii in sorted(set(control['impurity_problem_equivalence'])):
#         prob_ind=con3trol['impurity_problem_equivalence'].index(ii)
#         nimp_orb=len(control['impurity_wan'][prob_ind])    
#         transmat=np.zeros((nimp_orb,nimp_orb), dtype='complex')        
#         for jj in range(nimp_orb):
#             transmat2=array(map(float,g.readline().split()))
#             transmat[jj,:]=transmat2[0::2]+transmat2[1::2]*1j
#         trans_basis[str(ii)]=transmat

#     return trans_basis

# def read_impurity_vec_static(control,filename):
#     imp_basis={}
#     g=open(filename, 'r')

#     for ii in sorted(set(control['impurity_problem_equivalence'])):
#         prob_ind=control['impurity_problem_equivalence'].index(ii)
#         nimp_orb=len(control['impurity_wan'][prob_ind])        
#         impmat=np.zeros((nimp_orb,nimp_orb), dtype='complex')        
#         for jj in range(nimp_orb):
#             impmat2=array(map(float,g.readline().split()))
#             impmat[jj,:]=impmat2[0::2]+impmat2[1::2]*1j
#         imp_basis[str(ii)]=impmat
#     return imp_basis


def read_impurity_mat_static(control,filename):
    imp_basis={}
    g=open(filename, 'r')
    for ii in sorted(set(control['impurity_problem_equivalence'])):
        prob_ind=control['impurity_problem_equivalence'].index(ii)
        nimp_orb=len(control['impurity_wan'][prob_ind])

        impmat=np.zeros((nimp_orb,nimp_orb), dtype='complex')        
        # for jj in range(nimp_orb):
        #     impmat2=array([float(x) for x in g.readline().split()])
        #     for kk in range(0,nimp_orb*2,2):
        #         impmat[jj,kk]=impmat2[kk]+impmat2[kk+1]*1j
        for jj in range(nimp_orb):
            impmat2=array(list(map(float,g.readline().split())))
            impmat[jj,:]=impmat2[0::2]+impmat2[1::2]*1j
        imp_basis[str(ii)]=impmat
    return imp_basis


def read_impurity_mat_dynamic(control,filename):
    imp_basis={}
    dat=np.loadtxt(filename)
    print(np.shape(dat))

    start_array={}
    end_array={}

    last_index=1

    for ii in sorted(set(control['impurity_problem_equivalence'])):
        prob_ind=control['impurity_problem_equivalence'].index(ii)
        nimp_orb=len(control['impurity_wan'][prob_ind])
        start_array[ii]=last_index
        end_array[ii]=last_index+2*nimp_orb**2
        last_index=last_index+2*nimp_orb**2

    print(start_array)
    print(end_array)

    for ii in sorted(set(control['impurity_problem_equivalence'])):
        prob_ind=control['impurity_problem_equivalence'].index(ii)
        nimp_orb=len(control['impurity_wan'][prob_ind])
        dat3=np.reshape(dat[:,start_array[ii]:end_array[ii]], (control['n_omega'], 2, nimp_orb,nimp_orb), order='F')
        imp_basis[str(ii)]=dat3[:,0,:,:]+dat3[:,1,:,:]*1j

    return imp_basis



def cal_hyb_diagonal(control,imp):

    os.chdir(control['lowh_directory'])

    hyb_mat=read_impurity_mat_dynamic(control,control['lowh_directory']+'/delta_mat.dat')

    # print hyb_mat

    hyb_table=[]
    for jj in range(control['n_omega']):
        hyb_omega=[control['omega'][jj]]
        for ii in sorted(set(control['impurity_problem_equivalence'])):
            hyb_vec=imp_from_mat_to_array(hyb_mat[str(ii)][jj,:,:],imp[str(abs(ii))]['impurity_matrix'])

            hyb_omega=hyb_omega+np.reshape(np.stack((hyb_vec.real, hyb_vec.imag), 0), (len(hyb_vec)*2), order='F').tolist()
        hyb_table.append(hyb_omega)

    with open(control['lowh_directory']+'/delta.dat', 'w') as outputfile:
        outputfile.write(tabulate(hyb_table, headers=control['sig_header'], floatfmt=".12f", numalign="right",  tablefmt="plain"))

    if (control['method']=='lqsgw+dmft'):
        iter_string='_'+str(control['iter_num_impurity'])
    elif (control['method']=='lda+dmft'):
        iter_string='_'+str(control['iter_num_outer'])+'_'+str(control['iter_num_impurity'])
    labeling_file('./delta.dat', iter_string)

    shutil.copy('./delta.dat', control['top_dir'])    

    print('delta.dat generation done', file=control['h_log'],flush=True)
    causality=test_causality('./delta.dat')

    os.chdir(control['lowh_directory'])    

    return causality


# def cal_sig_dc_diagonal(control,imp):

#     os.chdir(control['dc_directory'])

#     trans_basis=read_impurity_mat_static(control,control['lowh_directory']+'/trans_basis.dat')
#     sig_mat=read_impurity_mat_dynamic(control,control['dc_directory']+'/delta_mat.dat')
#     h=open('./Delta.inp', 'w')

#     print hyb_mat

#     for jj in range(control['n_omega']):
#         h.write(str(control['omega'][jj])+'          ')
#         for ii in sorted(set(control['impurity_problem_equivalence'])):        
#             hyb_mat_new=dot(dot(trans_basis[str(ii)], hyb_mat[str(ii)][jj,:,:]), conj(transpose(trans_basis[str(ii)])))
#             hyb_vec=imp_from_mat_to_array(hyb_mat_new,imp[str(abs(ii))]['impurity_matrix'])
#             for kk in range(len(hyb_vec)):
#                 h.write(str(real(hyb_vec[kk]))+'   '+str(imag(hyb_vec[kk]))+'    ')
#         h.write('\n')
#     h.close()            
#     if (control['method']=='lqsgw+dmft'):
#         iter_string='_'+str(control['iter_num_impurity'])
#     elif (control['method']=='lda+dmft'):
#         iter_string='_'+str(control['iter_num_outer'])+'_'+str(control['iter_num_impurity'])
#     labeling_file('./Delta.inp', iter_string)

#     print('Delta.inp generation done', file=control['h_log'],flush=True)
#     causality=test_causality('./Delta.inp')

#     return causality


def labeling_file(filename,iter_string):
    dirname=os.path.abspath(os.path.dirname(filename))
    filenameonly=os.path.basename(filename)
    temp=filenameonly.split('.')
    shutil.copy(dirname+'/'+filenameonly, dirname+"/"+'.'.join(temp[0:-1])+iter_string+'.'+temp[-1])
    return None


def directory_setup(control):

    if (control['method'] =='lda+dmft'):
    #lattice
        tempdir=control['lattice_directory']
        if len(glob.glob(tempdir))==0 : os.mkdir(tempdir)
    elif  (control['method'] =='lqsgw+dmft'):
        tempdir=control['coulomb_directory']
        if len(glob.glob(tempdir))==0 : os.mkdir(tempdir)

    #wannier90 directory
    tempdir=control['wannier_directory']
    if len(glob.glob(tempdir))==0 : os.mkdir(tempdir)

    tempdir=control['dc_directory']
    if len(glob.glob(tempdir))==0 : os.mkdir(tempdir)

    # ctqmc
    tempdir=control['impurity_directory']
    if len(glob.glob(tempdir))==0 : os.mkdir(tempdir)    
    for ii in range(1,amax(control['impurity_problem_equivalence'])+1):
        tempdir=control['impurity_directory']+'/'+str(ii)
        if len(glob.glob(tempdir))==0 : os.mkdir(tempdir)        
        tempdir=control['dc_directory']+'/'+str(ii)
        if len(glob.glob(tempdir))==0 : os.mkdir(tempdir)        
    # delta
    tempdir=control['lowh_directory']
    if len(glob.glob(tempdir))==0 : os.mkdir(tempdir)

    return None

def check_for_files(filepath, h_log):
    if len(glob.glob(filepath))==0:
        print('missing:', filepath, file=control['h_log'],flush=True)
        quit()
    return None

def gaussian_broadening_linear(x, y, w1, temperature, cutoff):
    # broadening starts at the second matsubara points
    print(shape(x))
    print(shape(y))
    print(x)
    print(y)
    w0=(1.0-3.0*w1)*np.pi*temperature*8.6173303*10**-5
    width_array=w0+w1*x
    cnt=0
    ynew=np.zeros(len(y), dtype='complex')
    for x0 in x:
        if (x0>cutoff+(w0+w1*cutoff)*3.0):
            ynew[cnt]=y[cnt]
        else:
            if ((x0>3*width_array[cnt]) and ((x[-1]-x0)>3*width_array[cnt])):
                dist=1.0/np.sqrt(2*pi)/width_array[cnt]*exp(-(x-x0)**2/2.0/width_array[cnt]**2)
                ynew[cnt]=sum(dist*y)/sum(dist)
            else:
                ynew[cnt]=y[cnt]
        cnt=cnt+1
    return ynew

def solve_impurity_patrick(control):

    # execute CTQMC
    # chdir_string='cd '+control['top_dir']+'/impurity; '
    run_string=control['mpi_prefix_impurity']+' '+control['comsuitedir']+"/CTQMC params"
    cmd = run_string

    print(cmd, file=control['h_log'],flush=True)

    with open('./ctqmc.out', 'w') as logfile, open('./ctqmc.err', 'w') as errfile:
        ret = subprocess.call(cmd, shell=True,stdout = logfile, stderr = errfile)
        if ret != 0:
            print("Error in CTQMC. Check ctqmc.err for error message.", file=control['h_log'],flush=True)
            sys.exit()
    return None


def measure_impurity_patrick(control):

    measure_exe = control['comsuitedir'] + '/EVALSIM'

    run_string=measure_exe+' params'
    cmd = run_string

    print(cmd, file=control['h_log'],flush=True)

    with open('./evalsim.out', 'w') as logfile, open('./evalsim.err', 'w') as errfile :
        ret = subprocess.call(cmd,shell=True, stdout=logfile, stderr=errfile)
        if ret != 0:
            print("Error in EVALSIM. Check evalsim.err for error message.", file=control['h_log'],flush=True)
            sys.exit()
    print("measure self-energy done", file=control['h_log'],flush=True)
    if (control['method']=='lqsgw+dmft'):
        iter_string='_'+str(control['iter_num_impurity'])
    elif (control['method']=='lda+dmft'):
        iter_string='_'+str(control['iter_num_outer'])+'_'+str(control['iter_num_impurity'])    
    shutil.copy("./evalsim.out", "./evalsim"+iter_string+'.log')
    return None


def write_json_all(control,imp,data_array,json_name):
    # assume that it is diagonal matrix

    for key, value in imp.items(): # for the ordered phase this part should be fixed
        json_dict={}
        if (not (isinstance(imp[key], dict))):
            continue
        n_iio=amax(imp[key]['impurity_matrix'])
        if (imp[key]['para']):
            for kk in range(n_iio):
                orb_name=str(kk+1)
                json_dict[orb_name]={}
                json_dict[orb_name]['beta']=imp['beta']
                json_dict[orb_name]['real']=real(data_array[key][:,kk]).tolist()
                json_dict[orb_name]['imag']=imag(data_array[key][:,kk]).tolist()
        else:
            mkey=str(-int(key))
            for kk in range(n_iio):
                orb_name=str(kk+1)
                json_dict[orb_name]={}
                json_dict[orb_name]['beta']=imp['beta']
                json_dict[orb_name]['real']=real(data_array[key][:,kk]).tolist()
                json_dict[orb_name]['imag']=imag(data_array[key][:,kk]).tolist()

                orb_name=str(kk+1+n_iio)
                json_dict[orb_name]={}
                json_dict[orb_name]['beta']=imp['beta']
                json_dict[orb_name]['real']=real(data_array[mkey][:,kk]).tolist()
                json_dict[orb_name]['imag']=imag(data_array[mkey][:,kk]).tolist()

        with open(control['impurity_directory']+'/'+key+'/'+json_name,'w') as outfile:
            json.dump(json_dict, outfile,sort_keys=True, indent=4, separators=(',', ': '))
    print(json_name+" written", file=control['h_log'],flush=True)
    return None    



def read_json(jsonfile):
    Sig_temp=json.load(open(jsonfile))
    n_omega=len(Sig_temp['1']['real'])
    n_iio=len(Sig_temp.keys())
    dat1=zeros((n_omega, n_iio), dtype='complex')
    for key, value in Sig_temp.items():    
        dat1[:,int(key)-1]=array(Sig_temp[key]['real'])+array(Sig_temp[key]['imag'])*1j
    return dat1

def read_function_from_jsonfile(jsonfile, dict_name):
    Sig_temp=json.load(open(jsonfile))[dict_name]
    n_omega=len(Sig_temp['1']["function"]['real'])
    n_iio=len(Sig_temp.keys())
    dat1=zeros((n_omega, n_iio), dtype='complex')
    for key, value in Sig_temp.items():    
        dat1[:,int(key)-1]=array(Sig_temp[key]["function"]['real'])+array(Sig_temp[key]["function"]['imag'])*1j
    return dat1


def impurity_postprocessing(control, imp, key):

    if (control['method']=='lqsgw+dmft'):
        iter_string='_'+str(control['iter_num_impurity'])
    elif (control['method']=='lda+dmft'):
        iter_string='_'+str(control['iter_num_outer'])+'_'+str(control['iter_num_impurity'])

    labeling_file('./params.obs.json',iter_string)
    labeling_file('./params.meas.json',iter_string)    

    histo_temp=json.load(open('params.obs.json'))["expansion histogram"]

    histo=np.zeros((shape(histo_temp)[0], 2))
    histo[:,0]=arange(shape(histo_temp)[0])
    histo[:,1]=histo_temp
    nn=json.load(open('params.obs.json'))["scalar"]["N"]
    ctqmc_sign=json.load(open('params.obs.json'))["sign"]

    # histogram
    firstmoment=sum(histo[:,0]*histo[:,1])/sum(histo[:,1])
    secondmoment=sum((histo[:,0]-firstmoment)**2*histo[:,1])/sum(histo[:,1])
    thirdmoment=sum((histo[:,0]-firstmoment)**3*histo[:,1])/sum(histo[:,1])/secondmoment**(3.0/2.0)
    print('histogram information for impurity_'+imp['name'], file=control['h_log'],flush=True)
    print('first moment',  firstmoment,                      file=control['h_log'],flush=True)
    print('second moment', secondmoment,                     file=control['h_log'],flush=True)
    print('third moment',  thirdmoment,                      file=control['h_log'],flush=True)

    # previous_iter_string='_'.join(map(str,iter_string.split('_')[:-1]))+'_'+str(int(iter_string.split('_')[-1])-1)    


    green=read_function_from_jsonfile('./params.obs.json',"green")
    sigma_bare=read_function_from_jsonfile('./params.obs.json',"self-energy")

    sigma_old=array_impurity_dynamic(control,imp,control['impurity_directory']+'/sig.dat')
    sigma=np.zeros(np.shape(sigma_bare), dtype='complex')
    sigma_to_delta=np.zeros(np.shape(sigma_bare), dtype='complex')

    n_iio=amax(imp[key]['impurity_matrix'])

    sig_causality=1    

    for jj in range(n_iio):
        sigma[:,jj]=gaussian_broadening_linear(control['omega'], sigma_bare[:,jj], 0.05, imp['temperature'], imp[key]['green_cutoff'])
        if ((imag(sigma[:,jj])>0.0).any()):
            sig_causality=0
            sigma_to_delta[:,jj]=sigma_old[key][:,jj]
        else:
            sigma_to_delta[:,jj]=(sigma_old[key][:,jj])*(1.0-control['sigma_mix_ratio'])+(sigma[:,jj])*control['sigma_mix_ratio']
    if (not imp[key]['para']):
        for jj in range(n_iio, n_iio*2):
            mkey=str(-int(key))
            sigma[:,jj]=gaussian_broadening_linear(control['omega'], sigma_bare[:,jj], 0.05, imp['temperature'], imp[key]['green_cutoff'])
            if ((imag(sigma[:,jj])>0.0).any()):
                sig_causality=0
                sigma_to_delta[:,jj]=sigma_old[mkey][:,jj-n_iio]
            else:
                sigma_to_delta[:,jj]=(sigma_old[mkey][:,jj-n_iio])*(1.0-control['sigma_mix_ratio'])+(sigma[:,jj])*control['sigma_mix_ratio']

    if (imp[key]['para']):
        sig_diff_ave=np.sqrt(np.mean(np.absolute((sigma_to_delta-sigma_old[key]))**2))
    else:
        mkey=str(-int(key))
        sig_diff_ave=np.sqrt(np.mean((np.absolute((sigma_to_delta[:,0:n_iio]-sigma_old[key]))+np.absolute((sigma_to_delta[:,n_iio:]-sigma_old[mkey])))**2)/2.0)


    if (sig_causality==1):
        causality_flag='good'
    else:
        causality_flag='broken'

    if (control['method']=='lda+dmft'):
        control['conv_table'].append(['impurity_'+key,control['iter_num_outer'], '', control['iter_num_impurity'],causality_flag,'','','','',sig_diff_ave,nn,firstmoment,secondmoment,ctqmc_sign])
        with open(control['top_dir']+'/convergence.log', 'w') as outputfile:
            outputfile.write(tabulate(control['conv_table'], headers=control['convergence_header'], numalign="right",  floatfmt=".5f"))
    elif (control['method']=='lqsgw+dmft'):
        control['conv_table'].append(['impurity_'+key,control['iter_num_impurity'],causality_flag,'','','','',sig_diff_ave,nn,firstmoment,secondmoment,ctqmc_sign])
        with open(control['top_dir']+'/convergence.log', 'w') as outputfile:
            outputfile.write(tabulate(control['conv_table'], headers=control['convergence_header'], numalign="right",  floatfmt=".5f"))

    return green, sigma_bare, sigma, sigma_to_delta


def test_causality(filename):
    causality=1
    dat=np.loadtxt(filename)
    if ((dat[:,2::2]>0.0).any()):
        causality=0
        np.savetxt(filename+'b', dat)
        labeling_file(filename+'b',iter_string)
        print("Causality in "+filename+" is broken", file=control['h_log'],flush=True)
    else:
        print("Causality in "+filename+" is good",   file=control['h_log'],flush=True)
    return causality



def write_transformation_matrix(control, filename):
    os.chdir(control['lowh_directory'])
    

    if (control['trans_basis_mode']==2):
        f=open('trans_basis.dat', 'w')
        g=open(filename, 'r')
        
        for ii in sorted(set(control['impurity_problem_equivalence'])):
            prob_ind=control['impurity_problem_equivalence'].index(ii)
            nimp_orb=len(control['impurity_wan'][prob_ind])
            tempmat=np.zeros((nimp_orb,nimp_orb))    
            for jj in nimp_orb:
                tempmat[jj,:]=array(list(map(float,g.readline().split())))
            if (trace(tempmat) > control['metal_threshold']):
                w, v=np.linalg.eigh(tempmat)
                v=tranpose(v)
            else:
                v=np.identity(nimp_orb)
            for iorb in range(nimp_orb):    
                for jorb in range(nimp_orb):
                    f.write(str(v[iorb,jorb])+'     0.0     ')
                f.write("\n")

        f.close()
        g.close()
        shutil.copy('trans_basis.dat', control['top_dir'])

    if (control['method']=='lqsgw+dmft'):
        iter_string='_'+str(control['iter_num_impurity'])
    elif (control['method']=='lda+dmft'):
        iter_string='_'+str(control['iter_num_outer'])+'_'+str(control['iter_num_impurity'])        
    labeling_file('./trans_basis.dat', iter_string)

    os.chdir(control['top_dir'])    

    return None

def run_comlowh(control):

    os.chdir(control['lowh_directory'])
    run_string=control['mpi_prefix_lowh']+' '+control['comsuitedir']+"/ComLowH"
    logfilename=control['lowh_directory']+'/comlowh.out'
    errfilename=control['lowh_directory']+'/comlowh.err'    
    errormessage="Error in comlowh. Check comlowh.err for error message."

    cmd = run_string
    print(cmd, file=control['h_log'],flush=True)

    with open(logfilename, 'w') as logfile, open(errfilename, 'w') as errfile:
        ret = subprocess.call(cmd, shell=True,stdout = logfile, stderr = errfile)
        if ret != 0:
            print(errormessage, file=control['h_log'],flush=True)
            sys.exit()

    if (control['method']=='lqsgw+dmft'):
        iter_string="_"+str(control['iter_num_impurity'])
    elif (control['method']=='lda+dmft'):
        iter_string="_"+str(control['iter_num_outer'])+"_"+str(control['iter_num_impurity'])

    # labeling_file('./wannier_den_matrix.dat',iter_string)
    labeling_file('./comlowh.log',iter_string)
    labeling_file('./comlowh.out',iter_string)        
    labeling_file('./delta_mat.dat',iter_string)
    labeling_file('./g_loc_mat.dat',iter_string)
    labeling_file('./local_spectral_matrix_ef.dat',iter_string)
    labeling_file('./e_projected_mat.dat',iter_string)
    labeling_file('./ef.dat',iter_string)

    os.chdir(control['top_dir'])


    print("comlowh done", file=control['h_log'],flush=True)

    return None


def run_comcoulomb(control,imp):

    os.chdir(control['coulomb_directory'])
    run_string=control['mpi_prefix_coulomb']+' '+control['comsuitedir']+"/ComCoulomb"
    logfilename=control['coulomb_directory']+'/comcoulomb.out'
    errfilename=control['coulomb_directory']+'/comcoulomb.err'    
    errormessage="Error in comcomcoulomb. Check comcoulomb.err for error message."

    cmd = run_string
    print(cmd, file=control['h_log'],flush=True)

    with open(logfilename, 'w') as logfile, open(errfilename, 'w') as errfile:
        ret = subprocess.call(cmd, shell=True,stdout = logfile, stderr = errfile)
        if ret != 0:
            print(errormessage, file=control['h_log'],flush=True)
            sys.exit()

    iter_string="_"+str(control['iter_num_outer'])
    labeling_file('./comcoulomb.out',iter_string)
    labeling_file('./comcoulomb.ini',iter_string)
    files = glob.iglob(control['coulomb_directory']+"/*u_Slater*.rst")
    for filename in files:
        labeling_file(filename, iter_string)
    os.chdir(control['top_dir'])    

    return None


def comcoulomb_postprocessing(control,imp):
    slater_v={}
    slater_u={}
    slater_w={}    

    for ii in sorted(set(control['impurity_problem_equivalence'])):
        if (ii>0):
            jj=control['impurity_problem_equivalence'].index(ii)
            iatom=control['impurity_problem'][jj][0]
            shell=control['impurity_problem'][jj][1]

            if (shell=='s'):
                l_char='0'
            elif (shell=='p'):
                l_char='1'
            elif (shell=='d'):
                l_char='2'
            elif (shell=='f'):
                l_char='3'

            files = glob.iglob(control['coulomb_directory']+"/*_v_Slater_*"+str(iatom)+'_'+l_char+'.dat')
            for filename in files:
                slater_v[str(ii)]=np.loadtxt(filename)
                imp[str(ii)]['f0']=slater_v[str(ii)][0]
                if (int(l_char) >0):
                    imp[str(ii)]['f2']=slater_v[str(ii)][1]
                if (int(l_char) >1):    
                    imp[str(ii)]['f4']=slater_v[str(ii)][2]
                if (int(l_char) >2):    
                    imp[str(ii)]['f6']=slater_v[str(ii)][3]

            files = glob.iglob(control['coulomb_directory']+"/*_w_Slater_*"+str(iatom)+'_'+l_char+'.dat')
            for filename in files:
                tempmat=np.loadtxt(filename)

                n_nu=int(floor((tempmat[-1,0])/(2*pi/imp['beta'])))
                nu=arange(n_nu)*(2*pi/imp['beta'])


                dynamical_f0=cubic_interp1d(nu,tempmat[:,0], tempmat[:,1])
                if (int(l_char) >0):
                    dynamical_f2=cubic_interp1d(nu,tempmat[:,0], tempmat[:,2])
                if (int(l_char) >1):    
                    dynamical_f4=cubic_interp1d(nu,tempmat[:,0], tempmat[:,3])
                if (int(l_char) >2):    
                    dynamical_f6=cubic_interp1d(nu,tempmat[:,0], tempmat[:,4])

                if (int(l_char)==0):
                    slater_w[str(ii)]=transpose(vstack((dynamical_f0)))
                elif (int(l_char)==1):
                    slater_w[str(ii)]=transpose(vstack((dynamical_f0, dynamical_f2)))
                elif (int(l_char)==2):
                    slater_w[str(ii)]=transpose(vstack((dynamical_f0, dynamical_f2, dynamical_f4)))
                elif (int(l_char)==3):
                    slater_w[str(ii)]=transpose(vstack((dynamical_f0, dynamical_f2, dynamical_f4, dynamical_f6)))

            files = glob.iglob(control['coulomb_directory']+"/*_u_Slater_*"+str(iatom)+'_'+l_char+'.dat')
            for filename in files:

                tempmat=np.loadtxt(filename)

                n_nu=int(floor((tempmat[-1,0])/(2*pi/imp['beta'])))
                nu=arange(n_nu)*(2*pi/imp['beta'])

                dynamical_f0=cubic_interp1d(nu,tempmat[:,0], tempmat[:,1])
                if (int(l_char) >0):
                    dynamical_f2=cubic_interp1d(nu,tempmat[:,0], tempmat[:,2])
                if (int(l_char) >1):        
                    dynamical_f4=cubic_interp1d(nu,tempmat[:,0], tempmat[:,3])
                if (int(l_char) >2):        
                    dynamical_f6=cubic_interp1d(nu,tempmat[:,0], tempmat[:,4])    

                if (int(l_char)==0):
                    slater_u[str(ii)]=transpose(vstack((dynamical_f0)))
                elif (int(l_char)==1):
                    slater_u[str(ii)]=transpose(vstack((dynamical_f0, dynamical_f2)))
                elif (int(l_char)==2):
                    slater_u[str(ii)]=transpose(vstack((dynamical_f0, dynamical_f2, dynamical_f4)))
                elif (int(l_char)==3):
                    slater_u[str(ii)]=transpose(vstack((dynamical_f0, dynamical_f2, dynamical_f4, dynamical_f6)))

                imp[str(ii)]['dynamical_f0']=dynamical_f0-imp[str(ii)]['f0']

    u_table=nu
    w_table=nu
    # u_table=hstack((u_table, nu))
    # w_table=hstack((w_table, nu))    
    v_table=[]
    slater_header=['# nu(eV)']
    for ii in sorted(set(control['impurity_problem_equivalence'])):
        jj=control['impurity_problem_equivalence'].index(ii)
        iatom=control['impurity_problem'][jj][0]
        shell=control['impurity_problem'][jj][1]

        if (ii>0):
            if (shell=='s'):
                l_char='0'
            elif (shell=='p'):
                l_char='1'
            elif (shell=='d'):
                l_char='2'
            elif (shell=='f'):
                l_char='3'
                   
            u_table=column_stack((u_table, slater_u[str(ii)]))
            w_table=column_stack((w_table, slater_w[str(ii)]))
            v_table=hstack((v_table, slater_v[str(ii)]))

            slater_header.append(str(ii)+':f0(eV)')
            if (int(l_char)>0):                                            
                slater_header.append(str(ii)+':f2(eV)')
            if (int(l_char)>1):        
                slater_header.append(str(ii)+':f4(eV)')
            if (int(l_char)>2):            
                slater_header.append(str(ii)+':f6(eV)')
    with open(control['top_dir']+'/u_slater.dat', 'w') as outputfile:
        outputfile.write(tabulate(u_table, headers=slater_header, numalign="right",  floatfmt=".12f", tablefmt="plain"))
    with open(control['top_dir']+'/w_slater.dat', 'w') as outputfile:
        outputfile.write(tabulate(w_table, headers=slater_header, numalign="right",  floatfmt=".12f", tablefmt="plain"))
    slater_header=slater_header[1:]
    slater_header[0]='# '+slater_header[0]

    # print('v_table shape'+str(shape(v_table)), file=control['h_log'],flush=True)
    # print('v_table header shape'+str(shape(slater_header)), file=control['h_log'],flush=True)

    # print(v_table, file=control['h_log'],flush=True)
    # print(slater_header, file=control['h_log'],flush=True)    
    # print('v_table header shape'+str(shape(slater_header)), file=control['h_log'],flush=True)        
    
    with open(control['top_dir']+'/v_slater.dat', 'w') as outputfile:
        outputfile.write(tabulate([v_table], headers=slater_header, numalign="right",  floatfmt=".12f", tablefmt="plain"))

    print("comcoulomb done", file=control['h_log'],flush=True)

    return None



# def write_updates_json(control,imp):

#     if (control['spin_orbit']):
#         if (imp['problem']=='f'):
#             updates_json={
#                 "InsertEraseCSQ": {
#                 "Weight": 1.,
#                 "Moves": [
#                 [1.,"5/2,-5/2"],
#                 [1.,"5/2,-3/2"],
#                 [1.,"5/2,-1/2"],
#                 [1.,"5/2,+1/2"],
#                 [1.,"5/2,+3/2"],
#                 [1.,"5/2,+5/2"],
#                 [1.,"7/2,-7/2"],
#                 [1.,"7/2,-5/2"],
#                 [1.,"7/2,-3/2"],
#                 [1.,"7/2,-1/2"],
#                 [1.,"7/2,+1/2"],
#                 [1.,"7/2,+3/2"],
#                 [1.,"7/2,+5/2"],
#                 [1.,"7/2,+7/2"]        
#                 ]
#                 }
#                 }        
#     else:
#         if (imp['problem']=='d'):        
#             updates_json={
#                 "InsertEraseCSQ": {
#                 "Weight": 1.,
#                 "Moves": [
#                 [1., "yzUp"],
#                 [1., "zxUp"],
#                 [1., "xyUp"],
#                 [1., "3z2r2Up"],
#                 [1., "x2y2Up"],
#                 [1., "yzDown"],
#                 [1., "zxDown"],
#                 [1., "xyDown"],
#                 [1., "3z2r2Down"],
#                 [1., "x2y2Down"]
#                 ]
#                 }
#                 }
#     with open('Updates.json','w') as outfile:
#         json.dump(updates_json,outfile,sort_keys=True, indent=4, separators=(',', ': '))
#     print("Updates.json written"            , file=control['h_log'],flush=True)
#     return None

# def write_link_json(control, imp, key, equivalence_orb_mat):

#     # prob_ind=control['impurity_problem_equivalence'].index(int(key))
#     # nimp_orb=len(control['impurity_wan'][prob_ind])

#     if (control['spin_orbit']):
#         if (imp[key]['problem']=='f'):

#             link_json=[
#                 {
#                 "Irreps": ["5/2,-5/2"],
#                 "Flavors": [["5/2,-5/2"]],
#                 "Matrix": [
#                 ["+"+str(equivalence_orb_mat[0,0])+"+"]
#                 ]
#                 },
#                 {
#                 "Irreps": ["5/2,-3/2"],
#                 "Flavors": [["5/2,-3/2"]],
#                 "Matrix": [
#                 ["+"+str(equivalence_orb_mat[1,1])+"+"]                
#                 ]
#                 },
#                 {
#                 "Irreps": ["5/2,-1/2"],
#                 "Flavors": [["5/2,-1/2"]],
#                 "Matrix": [
#                 ["+"+str(equivalence_orb_mat[2,2])+"+"]                
#                 ]
#                 },
#                 {
#                 "Irreps": ["5/2,+1/2"],
#                 "Flavors": [["5/2,+1/2"]],
#                 "Matrix": [
#                 ["+"+str(equivalence_orb_mat[3,3])+"+"]                
#                 ]
#                 },
#                 {
#                 "Irreps": ["5/2,+3/2"],
#                 "Flavors": [["5/2,+3/2"]],
#                 "Matrix": [
#                 ["+"+str(equivalence_orb_mat[4,4])+"+"]                
#                 ]
#                 },
#                 {
#                 "Irreps": ["5/2,+5/2"],
#                 "Flavors": [["5/2,+5/2"]],
#                 "Matrix": [
#                 ["+"+str(equivalence_orb_mat[5,5])+"+"]                
#                 ]
#                 },
#                 {
#                 "Irreps": ["7/2,-7/2"],
#                 "Flavors": [["7/2,-7/2"]],
#                 "Matrix": [
#                 ["+"+str(equivalence_orb_mat[6,6])+"+"]                
#                 ]
#                 },
#                 {
#                 "Irreps": ["7/2,-5/2"],
#                 "Flavors": [["7/2,-5/2"]],
#                 "Matrix": [
#                 ["+"+str(equivalence_orb_mat[7,7])+"+"]                
#                 ]
#                 },
#                 {
#                 "Irreps": ["7/2,-3/2"],
#                 "Flavors": [["7/2,-3/2"]],
#                 "Matrix": [
#                 ["+"+str(equivalence_orb_mat[8,8])+"+"]                
#                 ]
#                 },
#                 {
#                 "Irreps": ["7/2,-1/2"],
#                 "Flavors": [["7/2,-1/2"]],
#                 "Matrix": [
#                 ["+"+str(equivalence_orb_mat[9,9])+"+"]                
#                 ]
#                 },
#                 {
#                 "Irreps": ["7/2,+1/2"],
#                 "Flavors": [["7/2,+1/2"]],
#                 "Matrix": [
#                 ["+"+str(equivalence_orb_mat[10,10])+"+"]                
#                 ]
#                 },
#                 {
#                 "Irreps": ["7/2,+3/2"],
#                 "Flavors": [["7/2,+3/2"]],
#                 "Matrix": [
#                 ["+"+str(equivalence_orb_mat[11,11])+"+"]                
#                 ]
#                 },
#                 {
#                 "Irreps": ["7/2,+5/2"],
#                 "Flavors": [["7/2,+5/2"]],
#                 "Matrix": [
#                 ["+"+str(equivalence_orb_mat[12,12])+"+"]                
#                 ]
#                 },
#                 {
#                 "Irreps": ["7/2,+7/2"],
#                 "Flavors": [["7/2,+7/2"]],
#                 "Matrix": [
#                 ["+"+str(equivalence_orb_mat[13,13])+"+"]                
#                 ]
#                 }
#                 ]
#     else:
#         if (imp[key]['problem']=='d'):
#             if (imp[key]['para']):
#                 index_shift=0
#             else:
#                 index_shift=amax(equivalence_orb_mat)            
#             link_json=[
#                 {
#                 "Irreps": ["yzUp"],
#                 "Flavors": [["yzUp"]],
#                 "Matrix": [
#                 ["+"+str(equivalence_orb_mat[0,0])+"+"]
#                 ]
#                 },
#                 {
#                 "Irreps": ["zxUp"],
#                 "Flavors": [["zxUp"]],
#                 "Matrix": [
#                 ["+"+str(equivalence_orb_mat[1,1])+"+"]                
#                 ]
#                 },
#                 {
#                 "Irreps": ["xyUp"],
#                 "Flavors": [["xyUp"]],
#                 "Matrix": [
#                 ["+"+str(equivalence_orb_mat[2,2])+"+"]                                
#                 ]
#                 },
#                 {
#                 "Irreps": ["3z2r2Up"],
#                 "Flavors": [["3z2r2Up"]],
#                 "Matrix": [
#                 ["+"+str(equivalence_orb_mat[3,3])+"+"]                                                
#                 ]
#                 },
#                 {
#                 "Irreps": ["x2y2Up"],
#                 "Flavors": [["x2y2Up"]],
#                 "Matrix": [
#                 ["+"+str(equivalence_orb_mat[4,4])+"+"]                                                                
#                 ]
#                 },
#                 {
#                 "Irreps": ["yzDown"],
#                 "Flavors": [["yzDown"]],
#                 "Matrix": [
#                 ["+"+str(equivalence_orb_mat[0,0]+index_shift)+"+"]                                                                
#                 ]
#                 },
#                 {
#                 "Irreps": ["zxDown"],
#                 "Flavors": [["zxDown"]],
#                 "Matrix": [
#                 ["+"+str(equivalence_orb_mat[1,1]+index_shift)+"+"]                                                                
#                 ]
#                 },
#                 {
#                 "Irreps": ["xyDown"],
#                 "Flavors": [["xyDown"]],
#                 "Matrix": [
#                 ["+"+str(equivalence_orb_mat[2,2]+index_shift)+"+"]                                                                
#                 ]
#                 },
#                 {
#                 "Irreps": ["3z2r2Down"],
#                 "Flavors": [["3z2r2Down"]],
#                 "Matrix": [
#                 ["+"+str(equivalence_orb_mat[3,3]+index_shift)+"+"]                                                                
#                 ]
#                 },
#                 {
#                 "Irreps": ["x2y2Down"],
#                 "Flavors": [["x2y2Down"]],
#                 "Matrix": [
#                 ["+"+str(equivalence_orb_mat[4,4]+index_shift)+"+"]                                                                
#                 ]
#                 }
#                 ]
#     with open('Link.json','w') as outfile:
#         json.dump(link_json,outfile,sort_keys=True, indent=4, separators=(',', ': '))
#     print("Link.json written"            , file=control['h_log'],flush=True)
#     return None


def write_params_json(control,imp,e_imp_key,trans_key,equivalence_key,beta):

    params_json={}
    # basis
    params_json["basis"]={}
    params_json["basis"]["orbitals"]=imp['problem'].lower()
    if (control['spin_orbit']):
        params_json["basis"]["type"]="coupled"
    else:
        params_json["basis"]["type"]="product"
    params_json["basis"]["transformation"]=trans_key.tolist()
    params_json["hloc"]={}

    mu_ctqmc=-e_imp_key[0,0]
    nimp_orb=len(imp['impurity_matrix'])
    e_ctqmc=(e_imp_key+identity(len(e_imp_key))*mu_ctqmc)

    params_json["hloc"]["one body"]=e_ctqmc.tolist()
    params_json["hloc"]["two body"]={}
    params_json["hloc"]["two body"]["parametrisation"]="slater-condon"    
    params_json["hloc"]["two body"]["F0"]=imp['f0']    
    if (params_json["basis"]["orbitals"]=='p') or (params_json["basis"]["orbitals"]=='d') or (params_json["basis"]["orbitals"]=='f') :
        params_json["hloc"]["two body"]["F2"]=imp['f2']
    if (params_json["basis"]["orbitals"]=='d') or (params_json["basis"]["orbitals"]=='f') :
        params_json["hloc"]["two body"]["F4"]=imp['f4']
    if (params_json["basis"]["orbitals"]=='f') :
        params_json["hloc"]["two body"]["F6"]=imp['f6']
        
    if imp["coulomb"]=="full":
        params_json["hloc"]["two body"]["approximation"]="none"
    elif imp["coulomb"]=="ising":
        params_json["hloc"]["two body"]["approximation"]="ising"
    params_json["hloc"]["quantum numbers"]={}
    params_json["hloc"]["quantum numbers"]["N"]={}    
    if (control['spin_orbit']):
        params_json["hloc"]["quantum numbers"]["Jz"]={}
    else:
        params_json["hloc"]["quantum numbers"]["Sz"]={}
    if (control['method']=='lqsgw+dmft'):
        params_json["dyn"]="dyn.json"
    params_json["beta"]=beta
    params_json["mu"]=mu_ctqmc
    params_json["hybridisation"]={}
    params_json["hybridisation"]["matrix"]=equivalence_key.tolist()
    params_json["hybridisation"]["functions"]="hyb.json"
    params_json["thermalisation time"]=imp['thermalization_time']
    params_json["measurement time"]=imp['measurement_time']
    params_json["green bulla"]=True    
    params_json["green matsubara cutoff"]=imp['green_cutoff']
    params_json["green basis"]="matsubara"

    # params_json["quantum number susceptibility"]=True
    # params_json["occupation susceptibility direct"]=True
    # params_json["occupation susceptibility bulla"]=True
    params_json["susceptibility cutoff"]=imp['susceptibility_cutoff']
    params_json["susceptibility tail"]=imp['susceptibility_tail']    
    # params_json["susceptibility tail"]=imp['susceptibility_cutoff']

    params_json["observables"]={}
    params_json["probabilities"]={}
    params_json["quantum numbers"]={}    
    if (control['spin_orbit']):
        params_json["observables"]["J2"]={}        
        params_json["probabilities"]=["N", "energy", "J2", "Jz"]
        params_json["quantum numbers"]["Jz"]={}        
    else:
        params_json["observables"]["S2"]={}                
        params_json["probabilities"]=["N", "energy", "S2", "Sz"]
        params_json["quantum numbers"]["Sz"]={}        

    params_json["occupation susceptibility bulla"]=True
    params_json["quantum number susceptibility"]=True

    params_json["occupation susceptibility direct"]=True
    params_json["occupation susceptibility bulla"]=True

    print("prepare_ctqmc:e_imp_done", file=control['h_log'],flush=True)

    with open('params.json','w') as outfile:
        json.dump(params_json,outfile, sort_keys=True, indent=4, separators=(',', ': '))
    print("params.json written", file=control['h_log'],flush=True)

    return None

def write_dynamical_f0_json(imp):

    with open('dyn.json','w') as outfile:
        json.dump(imp['dynamical_f0'].tolist(),outfile,sort_keys=True, indent=4, separators=(',', ': '))
    print("DynF0.json written"    , file=control['h_log'],flush=True)

    # os.chdir(control['top_dir'])
    return None



# def atom_run_patrick(control, imp):


#     # prob_ind=control['impurity_problem_equivalence'].index(int(key))
#     # nimp_orb=len(control['impurity_wan'][prob_ind])

#     if control['spin_orbit']:
#         if imp['problem']=='f':
#             atom_exe = control['comsuitedir'] + '/GA_F'        
#     else:
#         if imp['problem']=='d':
#             atom_exe = control['comsuitedir'] + '/GA_D'

#     # run_string=atom_exe+' params'
#     run_string='aprun -n 1 '+atom_exe+' params'    
#     cmd = run_string

#     print(cmd, file=control['h_log'],flush=True)

#     with open('./atom.out', 'w') as logfile:
#         ret = subprocess.call(cmd,shell=True, stdout=logfile, stderr=logfile)
#         if ret != 0:
#             print("Error in atom. Check atom.out for error message.", file=control['h_log'],flush=True)
#             sys.exit()            

#     print("prepare_ctqmc:atom done", file=control['h_log'],flush=True)
#     if (control['method']=='lqsgw+dmft'):
#         iter_string='_'+str(control['iter_num_impurity'])
#     elif (control['method']=='lda+dmft'):
#         iter_string='_'+str(control['iter_num_outer'])+'_'+str(control['iter_num_impurity'])    
#     shutil.copy("./atom.out", "./atom"+iter_string+'.log')
#     return None

def write_conv_dft(control):

    os.chdir(control['lattice_directory'])

    iter_string='_'+str(control['iter_num_outer'])
    f=open('./dft'+iter_string+'.out')
    cnt=0
    for line in f:
        temp=line.split()
        if (len(temp)==4):
            if temp[2]=='self-consistency=':
                cnt=cnt+1
                delta_rho=float(temp[3])
                control['conv_table'].append(['dft',control['iter_num_outer'],cnt,'', '', delta_rho, '','','','','','',''])
    with open(control['top_dir']+'/convergence.log', 'w') as outputfile:
        outputfile.write(tabulate(control['conv_table'], headers=control['convergence_header'], numalign="right",  floatfmt=".5f"))
    f.close()

    os.chdir(control['top_dir'])
    return None

def write_conv_coulomb(control,imp):

    os.chdir(control['coulomb_directory'])

    for ii in sorted(set(control['impurity_problem_equivalence'])):
        if (ii>0):
            control['conv_table'].append(['coulomb_'+str(ii),'', '', str(imp[str(ii)]['dynamical_f0'][0]+imp[str(ii)]['f0']), '','','','','','',''])
    with open(control['top_dir']+'/convergence.log', 'w') as outputfile:
        outputfile.write(tabulate(control['conv_table'], headers=control['convergence_header'], numalign="right",  floatfmt=".5f"))

    os.chdir(control['top_dir'])
    return None

def write_conv_wan(control):
    iter_string='_'+str(control['iter_num_outer'])    
    os.chdir(control['wannier_directory'])
    f=open('./wannier'+iter_string+'.wout')
    pp1=re.compile('Final State')
    cnt=0
    startline=0
    for line in f:
        mm1=pp1.search(line)
        if mm1:
            startline=cnt
        cnt=cnt+1 # start from 0
    f.close()

    f=open('./wannier'+iter_string+'.wout')
    lines=f.readlines()
    spmin=10000000.0
    spmax=0.0
    num_wann=np.shape(wan_hmat['basis'])[0]    
    wan_info=np.zeros((4,num_wann), order='F')
    cnt=0
    for ii in range(startline+1,startline+num_wann+1):
        wan_info[3,cnt]=float(lines[ii].split()[-1])
        temp1=lines[ii].split('(')[1]
        temp2=temp1.split(')')[0]
        # wan_info[:3,cnt]=[float(x) for x in temp2.split(',')]
        wan_info[:3,cnt]=list(map(float,temp2.split(',')))
        cnt=cnt+1
    f.close()

    # print wan_info

    f=open('./wannier'+iter_string+'.wout')
    lines=f.readlines()

    spmax=amax(wan_info[3,:])
    spmin=amin(wan_info[3,:])
    if (control['method']=='lda+dmft'):    
        control['conv_table'].append(['wannier',control['iter_num_outer'],'','','','', spmin,spmax,'','','','','',''])
        with open(control['top_dir']+'/convergence.log', 'w') as outputfile:
            outputfile.write(tabulate(control['conv_table'], headers=control['convergence_header'], numalign="right",  floatfmt=".5f"))
    if (control['method']=='lqsgw+dmft'):    
        control['conv_table'].append(['wannier','','','', spmin,spmax,'','','','','',''])
        with open(control['top_dir']+'/convergence.log', 'w') as outputfile:
            outputfile.write(tabulate(control['conv_table'], headers=control['convergence_header'], numalign="right",  floatfmt=".5f"))

    os.chdir(control['top_dir'])
    return None


def write_conv_delta(control,delta_causality):

    os.chdir(control['lowh_directory'])
    ef=np.float(np.loadtxt('ef.dat'))
    if (delta_causality==1):
        causality_flag='good'
    else:
        causality_flag='broken'

    if (control['method']=='lda+dmft'):
        control['conv_table'].append(['delta',control['iter_num_outer'],'',control['iter_num_impurity'],causality_flag,'','','', ef,'','','','',''])
        with open(control['top_dir']+'/convergence.log', 'w') as outputfile:
            outputfile.write(tabulate(control['conv_table'], headers=control['convergence_header'], numalign="right",  floatfmt=".5f"))
    if (control['method']=='lqsgw+dmft'):
        control['conv_table'].append(['delta',control['iter_num_impurity'],causality_flag,'','','', ef,'','','','',''])
        with open(control['top_dir']+'/convergence.log', 'w') as outputfile:
            outputfile.write(tabulate(control['conv_table'], headers=control['convergence_header'], numalign="right",  floatfmt=".5f"))
    os.chdir(control['top_dir'])
    return None


# def write_conv_imp(control,iter_string,iter_num_outer,iter_num_impurity,firstmoment,secondmoment,sig_causality,h_conv,h_log):

#     if (sig_causality==1):
#         causality_flag='good'
#     else:
#         causality_flag='broken'
#     os.chdir(control['impurity_directory'])
#     sig_ave=np.loadtxt('sig'+iter_string+'.dat')
#     sig=np.loadtxt('sig'+iter_string+'.dat')
#     sig_diff_ave=np.mean(np.absolute((sig_ave[:,1::2]+sig_ave[:,2::2]*1j)-(sig[:,1::2]+sig[:,2::2]*1j)))
#     nimp=read_nimp(imp_solver)
#     if (control['method']=='lda+dmft'):            
#         control['h_conv'].write('%1s%10s%10d%10s%10d%10s%10s%10s%10s%10s%10.7f%10.5f%10.3f%10.3f\n'%('','impurity',iter_num_outer,'',iter_num_impurity,causality_flag,'','','','',sig_diff_ave,nimp,firstmoment,secondmoment))
#     elif (control['method']=='lqsgw+dmft'):
#         control['h_conv'].write('%1s%10s%10d%10s%10s%10.7f%10.5f%10.3f%10.3f\n'%('','impurity',iter_num_impurity,causality_flag,'',sig_diff_ave,nimp,firstmoment,secondmoment))        
#     os.chdir(control['top_dir'])
#     return None                 

# def read_nimp(imp_solver):
#     # if imp_solver['solver']=='ctqmc_patrick':
#     nimp=np.loadtxt('N.dat')
#     # else:
#     #         f=open('sig.dat', 'r')
#     #         nimp=float((f.readline().split('=')[1]).split()[0])
#     #         f.close()
#     return nimp


def check_wannier_function_input(control,wan_hmat):    

    os.chdir(control['wannier_directory'])

    create_comwann_ini(control, wan_hmat)
    if os.path.exists(control['top_dir']+'/local_axis.dat'):
        shutil.copy(control['top_dir']+'/local_axis.dat', './')    

    os.chdir(control['top_dir'])
    return None


def check_coulomb_input(control):    

    os.chdir(control['coulomb_directory'])
    create_comcoulomb_ini(control)

    os.chdir(control['top_dir'])
    return None

def run_dft(control):

    os.chdir(control['lattice_directory'])
    iter_string='_'+str(control['iter_num_outer'])
    run_string=control['mpi_prefix_lattice']+' '+control['comsuitedir']+"/rspflapw.exe"
    cmd = run_string    

    with open(control['lattice_directory']+'/dft.out', 'w') as logfile, open(control['lattice_directory']+'/dft.err', 'w') as errfile:
        ret = subprocess.call(cmd, shell=True,stdout = logfile, stderr = errfile)
        if ret != 0:
            print("Error in dft. Check dft.err for error message.", file=control['h_log'],flush=True)
            sys.exit()

    allfile=control['allfile']

    labeling_file('./'+allfile+'.out',iter_string)    
    shutil.move('./dft.out', './dft'+iter_string+'.out')
    print("dft calculation done", file=control['h_log'],flush=True)
    os.chdir(control['top_dir'])
    return None

# def get_param_from_ini(param,stringstart,stringend,val_length,control):
#     f=open('ini', 'r')
#     pp=re.compile(param)
#     cnt=0
#     for line in f:
#         mm=pp.search(line)
#         if mm:
#             cnt=cnt+1
#             returnval=line[stringend:(stringend+val_length)]
#     if (cnt !=0):
#         return returnval.strip()
#     else:
#         print('couldn\'t find ', param, file=control['h_log'],flush=True)
#         quit()        


# def modify_chemical_potential_ubi(ef,h_log):
#     allfile=get_param_from_ini('allfile',1,10,72,control)            
#     allfile_out=string_addwhitespace(allfile, 72)
#     ef_old, ef_new=overwrite_rst.add_chemical_potential(allfile, 'dft', ef)
#     print('update, ef in dft', ef_old, ef_new, file=control['h_log'],flush=True)
#     return None

def prepare_dft_input(control):
    os.chdir(control['lattice_directory'])

    shutil.copy(control['lowh_directory']+"/wannier_den_matrix.dat", './')

    print("prepare_dft_input done", file=control['h_log'],flush=True)

    os.chdir(control['top_dir'])
    return None

# def  overwrite_restart_ubi(control):
#     f=open(control['allfile']+'.rst')
#     f.write('dft'+ '  0\n')
#     f.close()

# def check_nominal_dc_input(h_log):
#     check_for_files(control['top_dir']+'/dc/n_imp.dat', h_log)    

def cal_nominal_dc(imp,control):
    os.chdir(control['dc_directory'])
    f=open('dc_mat.dat', 'w')

    for ii in sorted(set(control['impurity_problem_equivalence'])):
        if (control['spin_orbit']):
            if (imp[str(abs(ii))]['problem']=='f'):
                nimp_orb=14
                uval=imp[str(abs(ii))]['f0']
                jval=(imp[str(abs(ii))]['f2']+imp[str(abs(ii))]['f4']+imp[str(abs(ii))]['f6'])/(6435.0/(286+195*0.668+250*0.494)*(1.0+0.668+0.494))
        else:
            if (imp[str(abs(ii))]['problem']=='f'):
                nimp_orb=7
                uval=imp[str(abs(ii))]['f0']
                jval=(imp[str(abs(ii))]['f2']+imp[str(abs(ii))]['f4']+imp[str(abs(ii))]['f6'])/(6435.0/(286+195*0.668+250*0.494)*(1.0+0.668+0.494))                
            elif (imp[str(abs(ii))]['problem']=='d'):
                nimp_orb=5
                uval=imp[str(abs(ii))]['f0']
                jval=(imp[str(abs(ii))]['f2']+imp[str(abs(ii))]['f4'])/14.0
            elif (imp[str(abs(ii))]['problem']=='p'):
                # from https://www.cond-mat.de/events/correl16/manuscripts/eder.pdf
                nimp_orb=3
                uval=imp[str(abs(ii))]['f0']
                jval=imp[str(abs(ii))]['f2']*5.0/25.0
            elif (imp[str(abs(ii))]['problem']=='s'):
                nimp_orb=1
                uval=imp[str(abs(ii))]['f0']
                jval=0.0
        dcval=(uval*(imp[str(abs(ii))]['nominal_n']-0.5)-jval*(imp[str(abs(ii))]['nominal_n']-1)*0.5)

        dcmat=np.identity(nimp_orb)*dcval
        for jj in range(nimp_orb):
            for kk in range(nimp_orb):
                f.write(str(dcmat[jj,kk])+'     0.0     ')
            f.write('\n')
    f.close()

    os.chdir(control['top_dir'])
    return None

def prepare_seed_dc_sig_and_wannier_dat(control,wan_hmat,imp):
    os.chdir(control['lowh_directory'])

    generate_comlowh_ini(control,wan_hmat,imp,1)

    natom=len(control['impurity_wan'])
    nimp_orb=0
    for ii in sorted(set(control['impurity_problem_equivalence'])):
        nimp_orb=nimp_orb+len(set(list(chain.from_iterable(imp[str(abs(ii))]['impurity_matrix'])))-{0})

    np.savetxt('dc.dat', zeros((1,nimp_orb*2)))

    aa=zeros((control['n_omega'],nimp_orb*2))
    bb=zeros((control['n_omega'],1))
    bb[:,0]=control['omega']
    np.savetxt('sig.dat',np.hstack((bb,aa)), header=' ')
    shutil.copy(control['wannier_directory']+"/wannier.dat", './')
    # make sig.dat

    os.chdir(control['top_dir'])
    return None

# def impurity_equivalence(control,imp):
#     imp_equivalence={}

#     num_atom=len(control['impurity_problem_equivalence'])
#     num_orb=zeros(num_atom, dtype=integer)
#     for ii in range(num_atom):
#         num_orb[ii]=len(control['impurity_wan'][ii])
#     iac=imp['impurity_atom_equivalence']
#     if (amin(iac) <0):
#         n_iac=amax(iac)*2
#         n_iac_nm=amax(iac)
#         n_iac_mat=n_iac+1
#         n_iac_mat_i=-n_iac_nm
#         n_iac_mat_f=n_iac_nm          
#         is_magnetic=1
#     else:
#         n_iac=amax(iac)
#         n_iac_nm=amax(iac)
#         n_iac_mat=n_iac
#         n_iac_mat_i=1
#         n_iac_mat_f=n_iac_nm                    
#         is_magnetic=0


#     num_orb_max=amax(num_orb)
#     ndeg_iac=zeros(n_iac_mat_f-n_iac_mat_i+1, dtype=integer)
#     norb_iac=zeros(n_iac_mat_f-n_iac_mat_i+1, dtype=integer)
#     ioac=zeros((num_orb_max,num_orb_max,n_iac_mat_f-n_iac_mat_i+1), dtype=integer)
#     n_ioac=amax(ioac)
#     iiiio=zeros((n_ioac,n_iac_mat_f-n_iac_mat_i+1), dtype=integer)
#     iio_diagonal=zeros((n_ioac,n_iac_mat_f-n_iac_mat_i+1), dtype=integer)
#     ndeg_ioac=zeros((n_ioac,n_iac_mat_f-n_iac_mat_i+1), dtype=integer)
#     ndeg_itot=zeros((n_ioac,n_iac_mat_f-n_iac_mat_i+1), dtype=integer)
#     ndeg_ioac_max=amax(ndeg_ioac)

#     for iatom in range(num_atom):
#         norb_iac[iac[iatom]-n_iac_mat_i]=num_orb[iatom]
#         ndeg_iac[iac[iatom]-n_iac_mat_i]=ndeg_iac[iac[iatom]-n_iac_mat_i]+1
#     for ii in (n_iac_mat_i, n_iac_mat_f):
#           if ((is_magnetic .eq. 1) .and. (ii .eq. 0)) cycle
#           do iorb=1, norb_iac(ii)
#             read(10,*) (ioac(iorb,jorb,ii),
#      $        jorb=1, norb_iac(ii))
#           enddo
#         enddo


def generate_comlowh_ini(control,wan_hmat,imp,is_recal_ef):

    f=open('comlowh.ini', 'w')
    f.write('1\n')
    natom=len(control['impurity_wan'])
    # nimp_orb=np.shape(control['impurity_wan'])[1]
    nimp_orb=np.zeros(natom, dtype=integer)
    for ii in range(natom):
        nimp_orb[ii]=len(control['impurity_wan'][ii])
    f.write(str(natom)+'\n')
    f.write(' '.join(map(str,nimp_orb))+'\n')
    f.write(' '.join(map(str,control['impurity_problem_equivalence']))+'\n')
    for ii in sorted(set(control['impurity_problem_equivalence'])):
        prob_ind=control['impurity_problem_equivalence'].index(ii)
        nimp_orb=len(control['impurity_wan'][prob_ind])        
        for jj in range(nimp_orb):
            f.write(' '.join(map(str,imp[str(abs(ii))]['impurity_matrix'][jj]))+'\n')

    for iatom in range(natom):
        f.write(' '.join(map(str,control['impurity_wan'][iatom]))+' ')
    f.write('\n')
    f.write(str(control['proj_win_min'])+'   '+str(control['proj_win_max'])+'\n')
    n_omega=control['n_omega']
    f.write(str(n_omega)+'\n')
    f.write('0.0\n')
    f.write('0.0\n')
    f.write(str(imp['beta'])+'\n')
    f.write(str(control['doping'])+'\n')
    f.write(str(is_recal_ef)+'\n')    
    f.write('bnd\n')    
    if (control['spin_orbit']):
        f.write('1\n')
    else:
        f.write('0\n')
    # if (control['update_mu_dmft_scf']):
    #         f.write('1\n')
    # else:
    #         f.write('0\n')
    f.write(' '.join(map(str,wan_hmat['kgrid']))+'\n')

    f.close()
    return None


def prepare_dc(control,wan_hmat,imp):
    if ('dc_mat_to_read' not in control):
        if (control['method']=='lqsgw+dmft'):
            gloc_mat=read_impurity_mat_dynamic(control,control['lowh_directory']+'/g_loc_mat.dat')
            trans_basis=read_impurity_mat_static(control,control['lowh_directory']+'/trans_basis.dat')
            print(trans_basis)
            for key, value in imp.items(): # for the ordered phase this part should be fixed
                if (not (isinstance(imp[key], dict))):
                    continue
                nimp_orb=len(imp[key]['impurity_matrix'])
                os.chdir(control['dc_directory']+'/'+key)
                f=open('comdc.ini', 'w')
                f.write(str(nimp_orb)+'\n')
                if (control['spin_orbit']):
                    f.write('1\n')
                else:
                    f.write('0\n')
                    f.write('0\n')
                f.close()
                f=open('g_loc.dat', 'w')
                for ii in range(control['n_omega']):
                    f.write(str(control['omega'][ii])+'  '+' '.join(map("{:.12f}".format, reshape(np.stack((gloc_mat[key][ii,:,:].real,gloc_mat[key][ii,:,:].imag),0), (2*nimp_orb**2), order='F')))+'\n')
                f.close()
                savetxt('trans_dc.dat',reshape(np.stack((trans_basis[key].real,trans_basis[key].imag),-1), (nimp_orb, 2*nimp_orb)))
                f=open('slater.dat', 'w')
                if (imp[key]['problem']=='s'):
                    f.write(str(imp[key]['f0'])+'\n')
                elif (imp[key]['problem']=='p'):
                    f.write(str(imp[key]['f0'])+'  '+str(imp[key]['f2'])+'\n')    
                elif (imp[key]['problem']=='d'):
                    f.write(str(imp[key]['f0'])+'  '+str(imp[key]['f2'])+'  '+str(imp[key]['f4'])+'\n')
                elif (imp[key]['problem']=='f'):
                    f.write(str(imp[key]['f0'])+'  '+str(imp[key]['f2'])+'  '+str(imp[key]['f4'])+'  '+str(imp[key]['f6'])+'\n')
                f.close()
                for ii in range(len(imp[str(key)]['dynamical_f0'])):
                    if (abs(imp[str(key)]['dynamical_f0'][ii]) <0.5):
                        break
                savetxt('dynamical_f0.dat', imp[str(key)]['dynamical_f0'][:ii])
            os.chdir(control['top_dir'])
    return None

def write_conv_dc(control,imp):

    if (control['method']=='lqsgw+dmft'):
        for key, value in imp.items(): # for the ordered phase this part should be fixed
            if (not (isinstance(imp[key], dict))):
                continue
            os.chdir(control['dc_directory']+'/'+key)

            # nimp_orb=len(imp[key]['impurity_matrix'])
            # nimp=np.trace(np.reshape(np.loadtxt('nimp.dat'), (nimp_orb,nimp_orb,5), order='F')[:,:,3])
            control['conv_table'].append(['dc_'+key,'','good','','','','','','','','',''])
            with open(control['top_dir']+'/convergence.log', 'w') as outputfile:
                outputfile.write(tabulate(control['conv_table'], headers=control['convergence_header'], numalign="right",  floatfmt=".5f"))
    return None

def run_dc(control,imp):
    if ('dc_mat_to_read' in control):
        os.chdir(control['dc_directory'])
        shutil.copy(control['dc_mat_to_read'], './dc_mat.dat')
        os.chdir(control['top_dir'])
    else:
        if (control['method']=='lda+dmft'):
            cal_nominal_dc(imp,control)
        elif (control['method']=='lqsgw+dmft'):
            for key, value in imp.items(): # for the ordered phase this part should be fixed
                if (not (isinstance(imp[key], dict))):
                    continue
                os.chdir(control['dc_directory']+'/'+key)

                run_string=control['mpi_prefix_dc']+' '+control['comsuitedir']+"/ComDC"
                logfilename=control['dc_directory']+'/'+key+'/comdc.out'
                errfilename=control['dc_directory']+'/'+key+'/comdc.err'
                errormessage="Error in comdc. Check comdc.err for error message."

                cmd = run_string
                print(cmd, file=control['h_log'],flush=True)

                with open(logfilename, 'w') as logfile, open(errfilename, 'w') as errfile:
                    ret = subprocess.call(cmd, shell=True,stdout = logfile, stderr = errfile)
                    if ret != 0:
                        print(errormessage, file=control['h_log'],flush=True)
                        sys.exit()
                iter_string="_"+str(control['iter_num_outer'])
                labeling_file('./comdc.out',iter_string)
                labeling_file('./sig_mat.dat',iter_string)

            os.chdir(control['dc_directory'])
            f=open(control['dc_directory']+'/dc_mat.dat','w')
            g=open(control['dc_directory']+'/zinv_m1_mat.dat','w')
            for ii in sorted(set(control['impurity_problem_equivalence'])):
                nimp_orb=len(imp[str(abs(ii))]['impurity_matrix'])
                dc=reshape(loadtxt(control['dc_directory']+'/'+str(abs(ii))+'/sig_mat.dat')[0,1:], (2,nimp_orb,nimp_orb), order='F')
                for jj in range(nimp_orb):
                    for kk in range(nimp_orb):
                        f.write(str(dc[0,jj,kk])+'     0.0     ')
                        g.write(str(-dc[1,jj,kk]/control['omega'][0])+'     0.0     ')
                    f.write('\n')
                    g.write('\n')
            f.close()
            g.close()

            sig_dc={}

            for ii in sorted(set(control['impurity_problem_equivalence'])):
                nimp_orb=len(imp[str(abs(ii))]['impurity_matrix'])                
                tempdat=reshape(loadtxt(control['dc_directory']+'/'+str(abs(ii))+'/sig_mat.dat')[:,1:], (control['n_omega'],2,nimp_orb,nimp_orb), order='F')
                sig_dc[str(ii)]=tempdat[:,0,:,:]+tempdat[:,1,:,:]*1j

            sig_table=[]
            for jj in range(control['n_omega']):
                sig_omega=[control['omega'][jj]]
                for ii in sorted(set(control['impurity_problem_equivalence'])):
                    sig_dc_vec=imp_from_mat_to_array(sig_dc[str(ii)][jj,:,:],imp[str(abs(ii))]['impurity_matrix'])
                    print(sig_dc_vec, np.reshape(np.stack((sig_dc_vec.real, sig_dc_vec.imag), 0), (len(sig_dc_vec)*2), order='F').tolist())
                    sig_omega=sig_omega+np.reshape(np.stack((sig_dc_vec.real, sig_dc_vec.imag), 0), (len(sig_dc_vec)*2), order='F').tolist()
                sig_table.append(sig_omega)


            with open(control['top_dir']+'/sig_dc.dat', 'w') as outputfile:
                outputfile.write(tabulate(sig_table, headers=control['sig_header'], floatfmt=".12f", numalign="right",  tablefmt="plain"))

            sig_hf_dc={}

            for ii in sorted(set(control['impurity_problem_equivalence'])):
                nimp_orb=len(imp[str(abs(ii))]['impurity_matrix'])                
                tempdat=reshape(loadtxt(control['dc_directory']+'/'+str(abs(ii))+'/hartree.dat')[:,3:], (nimp_orb,nimp_orb,2), order='F')+reshape(loadtxt(control['dc_directory']+'/'+str(abs(ii))+'/exchange.dat')[:,3:], (nimp_orb,nimp_orb,2), order='F')
                sig_hf_dc[str(ii)]=tempdat[:,:,0]+tempdat[:,:,1]*1j

            sig_table=[]
            hf_header=control['sig_header'][1:]
            hf_header[0]='# '+hf_header[0]
            for ii in sorted(set(control['impurity_problem_equivalence'])):
                dc_vec=imp_from_mat_to_array(sig_hf_dc[str(ii)],imp[str(abs(ii))]['impurity_matrix'])
                sig_table.append(np.reshape(np.stack((dc_vec.real, dc_vec.imag), 0), (len(dc_vec)*2), order='F').tolist())

            with open(control['top_dir']+'/sig_dc_hf.dat', 'w') as outputfile:
                outputfile.write(tabulate(sig_table, headers=hf_header, floatfmt=".12f", numalign="right",  tablefmt="plain"))


            labeling_file('./dc_mat.dat',iter_string)
            labeling_file('./zinv_m1_mat.dat',iter_string)
            os.chdir(control['top_dir'])
    return None

def generate_initial_transformation(control):
    os.chdir(control['lowh_directory'])
    print(control['impurity_wan'], file=control['h_log'],flush=True)

    if (control['trans_basis_mode']==0):
        f=open('trans_basis.dat', 'w')
        for ii in sorted(set(control['impurity_problem_equivalence'])):
            prob_ind=control['impurity_problem_equivalence'].index(ii)
            nimp_orb=len(control['impurity_wan'][prob_ind])
            transmat=np.identity(nimp_orb)
            for jj in range(nimp_orb):    
                for kk in range(nimp_orb):
                    f.write(str(transmat[jj,kk])+'     0.0      ')
                f.write("\n")
        f.close()
    elif (control['trans_basis_mode']==1):
        shutil.copy(control['trans_basis'], './trans_basis.dat')
    elif (control['trans_basis_mode']==2):
        if ('trans_basis' in control):
            shutil.copy(control['trans_basis'], './trans_basis.dat')
        else:
            f=open('trans_basis.dat', 'w')
            for ii in sorted(set(control['impurity_problem_equivalence'])):
                prob_ind=control['impurity_problem_equivalence'].index(ii)
                nimp_orb=len(control['impurity_wan'][prob_ind])
                transmat=np.identity(nimp_orb)
                for jj in range(nimp_orb):    
                    for kk in range(nimp_orb):
                        f.write(str(transmat[jj,kk])+'     0.0      ')
                    f.write("\n")
            f.close()

    shutil.copy('trans_basis.dat', control['top_dir'])
    
    if (control['method']=='lqsgw+dmft'):
        iter_string='_0'
    elif (control['method']=='lda+dmft'):
        iter_string='_1_0'
    labeling_file('./trans_basis.dat',iter_string)
    os.chdir(control['top_dir'])
    return None


def prepare_comlowh(control,wan_hmat,imp):

    os.chdir(control['lowh_directory'])

    generate_comlowh_ini(control,wan_hmat,imp,1)

    # wannier files
    shutil.copy(control['wannier_directory']+"/wannier.dat", './')
    shutil.copy(control['dc_directory']+"/dc.dat", './')
    shutil.copy(control['impurity_directory']+"/sig.dat", './')
    # if (control['dc_linear']):
    files = glob.iglob(control['dc_directory']+"/zinv_m1.dat")
    for filename in files:
        shutil.copy(filename, './')



    print("check and preparation done for the calculation of delta", file=control['h_log'],flush=True)
    os.chdir(control['top_dir'])
    return None

def comwann_postprocessing(control, wan_hmat):
    wan_hmat['basis']=read_wan_hmat_basis(control)
    find_impurity_wan(control, wan_hmat)    
    return None

def run_comwann(control,wan_hmat):

    os.chdir(control['wannier_directory'])

    run_string=control['mpi_prefix_wannier']+' '+control['comsuitedir']+"/ComWann"
    cmd = run_string
    print(cmd, file=control['h_log'],flush=True)

    with open(control['wannier_directory']+'/comwann.out', 'w') as logfile, open(control['wannier_directory']+'/comwann.err', 'w') as errfile:
        ret = subprocess.call(cmd,shell=True, stdout = logfile, stderr = errfile)
        if ret != 0:
            print("Error in comwann. Check comwann.err or OUT.", file=control['h_log'],flush=True)
            sys.exit()
    # shutil.move('./wannier_1.wout','./wannier.wout')    

    iter_string='_'+str(control['iter_num_outer'])

    labeling_file('./wannier.dat',iter_string)
    labeling_file('./wannier.chk',iter_string)
    labeling_file('./wannier.inip',iter_string)
    labeling_file('./wannier.eig',iter_string)
    labeling_file('./wannier.win',iter_string)
    labeling_file('./orb_for_froz_win.dat',iter_string)                
    shutil.copy('./wannier.wout','./wannier'+iter_string+'.wout')

    os.chdir(control['top_dir'])    

    return None

def cubic_interp1d(x0, x, y):
    """
    Interpolate a 1-D function using cubic splines.
        x0 : a float or an 1d-array
        x : (N,) array_like
            A 1-D array of real/complex values.
        y : (N,) array_like
            A 1-D array of real values. The length of y along the
            interpolation axis must be equal to the length of x.

    Implement a trick to generate at first step the cholesky matrice L of
    the tridiagonal matrice A (thus L is a bidiagonal matrice that
    can be solved in two distinct loops).

    additional ref: www.math.uh.edu/~jingqiu/math4364/spline.pdf 
    """
    x = np.asfarray(x)
    y = np.asfarray(y)

    # remove non finite values
    # indexes = np.isfinite(x)
    # x = x[indexes]
    # y = y[indexes]

    # check if sorted
    if np.any(np.diff(x) < 0):
        indexes = np.argsort(x)
        x = x[indexes]
        y = y[indexes]

    size = len(x)

    xdiff = np.diff(x)
    ydiff = np.diff(y)

    # allocate buffer matrices
    Li = np.empty(size)
    Li_1 = np.empty(size-1)
    z = np.empty(size)

    # fill diagonals Li and Li-1 and solve [L][y] = [B]
    Li[0] = sqrt(2*xdiff[0])
    Li_1[0] = 0.0
    B0 = 0.0 # natural boundary
    z[0] = B0 / Li[0]

    for i in range(1, size-1, 1):
        Li_1[i] = xdiff[i-1] / Li[i-1]
        Li[i] = sqrt(2*(xdiff[i-1]+xdiff[i]) - Li_1[i-1] * Li_1[i-1])
        Bi = 6*(ydiff[i]/xdiff[i] - ydiff[i-1]/xdiff[i-1])
        z[i] = (Bi - Li_1[i-1]*z[i-1])/Li[i]

    i = size - 1
    Li_1[i-1] = xdiff[-1] / Li[i-1]
    Li[i] = sqrt(2*xdiff[-1] - Li_1[i-1] * Li_1[i-1])
    Bi = 0.0 # natural boundary
    z[i] = (Bi - Li_1[i-1]*z[i-1])/Li[i]

    # solve [L.T][x] = [y]
    i = size-1
    z[i] = z[i] / Li[i]
    for i in range(size-2, -1, -1):
        z[i] = (z[i] - Li_1[i-1]*z[i+1])/Li[i]

    # find index
    index = x.searchsorted(x0)
    np.clip(index, 1, size-1, index)

    xi1, xi0 = x[index], x[index-1]
    yi1, yi0 = y[index], y[index-1]
    zi1, zi0 = z[index], z[index-1]
    hi1 = xi1 - xi0

    # calculate cubic
    f0 = zi0/(6*hi1)*(xi1-x0)**3 + \
         zi1/(6*hi1)*(x0-xi0)**3 + \
         (yi1/hi1 - zi1*hi1/6)*(x0-xi0) + \
         (yi0/hi1 - zi0*hi1/6)*(xi1-x0)
    return f0


# def modify_ini(flag,stringstart,stringend,val_length,val):
#     f=open('ini','r')
#     g=open('ini_new','w')
#     # print flag,stringstart,stringend,val_length,val
#     pp=re.compile(flag)
#     len_val=len(val)
#     len_flag=len(flag)

#     # print len_val, len_flag

#     stringlength=stringend-stringstart+1
#     newflag=flag+' '*(stringlength-len_flag-1)+'='+' '*(val_length-len_val)+val
#     # print stringlength
#     # print newflag

#     for line in f:
#         mm=pp.search(line)        
#         if mm:
#             # print line[:(stringstart-1)]
#             # print line[(stringend+val_length-1):]
#             newline=line[:(stringstart-1)]+newflag+line[(stringend+val_length):].rstrip()+'\n'
#             # print newline
#             g.write(newline)
#         else:
#             g.write(line)
#     g.close()
#     f.close()
#     shutil.move("ini_new", "ini")
#     return None



def optimized_nproc_for_comcoulomb(var1,npnt,ntau,nomega,nnu):
#     mpicom,flag,snproc=var1.split(" ")
    svar={}
    svar=var1.split(" ")
#     print "svar",svar
    snproc=""
    if var1.find("srun") != -1 :
        mpicom="srun"
        if var1.find("-np") != -1 :
            flag="-np"
            count=0
            for x in svar :
                if x==flag :
                    snproc=svar[count+1]
                    break
                count=count+1
            if snproc == "" :
                print("Error on finding nproce")
                exit()
        elif var1.find("-n") != -1 :
            flag="-n"
            count=0
            for x in svar :
                if x==flag :
                    snproc=svar[count+1]
                    break
                count=count+1
            if snproc == "" :
                print("Error on finding nproce")
                exit()         
        else :
            print("Error on finding -n or -np")
            exit()  

    elif var1.find("aprun") != -1 :
        mpicom="aprun"
        if var1.find("-np") != -1 :
            flag="-np"
            count=0
            for x in svar :
                if x==flag :
                    snproc=svar[count+1]
                    break
                count=count+1
            if snproc == "" :
                print("Error on finding nproce")
                exit()
        elif var1.find("-n") != -1 :
            flag="-n"
            count=0
            for x in svar :
                if x==flag :
                    snproc=svar[count+1]
                    break
                count=count+1
            if snproc == "" :
                print("Error on finding nproce")
                exit()
        else :
            print("Error on finding -n or -np")
            exit()

    elif  var1.find("mpirun") != -1 :
        mpicom="mpirun"
        if var1.find("-np") != -1 :
            flag="-np"
            count=0
            for x in svar :
                if x==flag :
                    snproc=svar[count+1]
                    break
                count=count+1
            if snproc == "" :
                print("Error on finding nproce")
                exit()
        elif var1.find("-n") != -1 :
            flag="-n"
            count=0
            for x in svar :
                if x==flag :
                    snproc=svar[count+1]
                    break
                count=count+1
            if snproc == "" :
                print("Error on finding nproce")
                exit()
        else :
            print("Error on finding -n or -np")
            exit()

    elif  var1.find("mpiexec") != -1 :
        mpicom="mpiexec"
        if var1.find("-np") != -1 :
            flag="-np"
            count=0
            for x in svar :
                if x==flag :
                    snproc=svar[count+1]
                    break
                count=count+1
            if snproc == "" :
                print("Error on finding nproce")
                exit()
        elif var1.find("-n") != -1 :
            flag="-n"
            count=0
            for x in svar :
                if x==flag :
                    snproc=svar[count+1]
                    break
                count=count+1
            if snproc == "" :
                print("Error on finding nproce")
                exit()
        else :
            print("Error on finding -n or -np")
            exit()

    else :
        print("Error on finding mpi command")
        exit()
#     snproc="400"

    nproc=int(snproc)
    Ntau=int(ntau/2+1)
    Nomega=nomega+1
    Nnu=nnu+1
    Nsmallest=Ntau
    if Nomega < Nsmallest :
        Nsmallest=Nomega
    if Nnu < Nsmallest :
        Nsmallest=Nnu
#     print 'Nsmallest',Nsmallest,'Ntau',Ntau,'Nomega',Nomega,'Nnu',Nnu
    NCom=0
    ListCom={}
    for i in range(1,Nsmallest+1):
        ltf=True
        if Ntau%i != 0 :
            ltf=False
        if Nomega%i != 0 :
            ltf=False
        if Nnu%i != 0 :
            ltf=False
        if ltf == True :
            NCom = NCom + 1
            ListCom[NCom]=i

#     print 'NCom',NCom,'ListCom',ListCom
    Sol=[[0 for x in range(3)] for y in range(NCom)]
    for i in range(0,NCom):
        ValCom=ListCom[i+1]
        ValQ=nproc//ValCom
        Sol[i][0]=ValCom
        Sol[i][1]=ValQ
        if ValQ < npnt :
            Sol[i][2]=nproc-ValCom*ValQ
        else:
            Sol[i][2]=nproc-ValCom*npnt
#     for i in range(0,NCom) :
#       print 'sol=',i,'nproc_tau',Sol[i][0],'nproc_k',Sol[i][1],'waste',Sol[i][2]
#     print ""
    Sol.sort(key=lambda x: x[2])    
#     for i in range(0,NCom) :
#       print 'sol=',i,'nproc_tau',Sol[i][0],'nproc_k',Sol[i][1],'waste',Sol[i][2]

    startval=Sol[0][2]
    istart=0
    iend=0
    for i in range(1,NCom):
        if Sol[i][2] != startval :
            iend=i-1
            break

#     print ""
    SolPart=[[0 for x in range(3)] for y in range(iend+1)]
    for i in range(0,iend+1):
        SolPart[i]=Sol[i]
    SolPart.sort(key=lambda x: x[0], reverse=True)
#     for i in range(0,iend+1) :
#       print 'sol=',i,'nproc_tau',SolPart[i][0],'nproc_k',SolPart[i][1],'waste',SolPart[i][2]


    nproc_tau=SolPart[0][0]
    nproc_k=SolPart[0][1]
    nproc0=nproc_tau*nproc_k

    varcom=mpicom+" "+flag+" "+str(nproc0)

    return varcom,nproc_k,nproc_tau


def find_allfile(dft_dir):
    files = glob.iglob(dft_dir+"/*.rst")
    allfile = None
    for filename in files:
        temp=filename[:-4].split('/')[-1].split('_')
        if len(temp)==1:
            if (temp[0]!='info' and temp[0]!='info2'):
                allfile=temp[0]
    return allfile


def lda_dmft(control,wan_hmat,imp):

    print("\n", file=control['h_log'],flush=True)
    print("\n", file=control['h_log'],flush=True)
    print("\n", file=control['h_log'],flush=True)

    while control['iter_num_outer'] <= control['max_iter_num_outer']:

        # iter_string_outer="_"+str(iter_num_outer)

        print("************************************************", file=control['h_log'],flush=True)
        print("iteration:  "+str(control['iter_num_outer']),      file=control['h_log'],flush=True)
        print("************************************************", file=control['h_log'],flush=True)

        control['iter_num_impurity']=0

        if (control['iter_num_outer']==1):
            initial_lattice_directory_setup(control)
        else:
            prepare_dft_input(control)
            run_dft(control)
            write_conv_dft(control)

        print("wannier function construction", file=control['h_log'],flush=True)
        if control['iter_num_outer']==1:
            prepare_initial_ef(control)

        check_wannier_function_input(control,wan_hmat)
        run_comwann(control, wan_hmat)
        comwann_postprocessing(control, wan_hmat)
        write_conv_wan(control)

        if control['iter_num_outer']==1:
            generate_initial_transformation(control)
            run_dc(control,imp)
            generate_initial_transformation(control)
            cal_dc_diagonal(control)
            generate_initial_self_energy(control,imp)



        control['iter_num_impurity']=1

        while control['iter_num_impurity'] <= control['max_iter_num_impurity']:

            print("\n",                                                                        file=control['h_log'],flush=True)
            print('*****   iter_num_impurity: ', str(control['iter_num_impurity']), '  *****', file=control['h_log'],flush=True)

            prepare_comlowh(control,wan_hmat,imp)
            run_comlowh(control)

            delta_causality=delta_postprocessing(control,imp)

            write_conv_delta(control,delta_causality)


            prepare_impurity_solver(control,wan_hmat,imp)
            run_impurity_solver(control,imp)

            control['iter_num_impurity']=control['iter_num_impurity']+1

        print("\n", file=control['h_log'],flush=True)
        print("\n", file=control['h_log'],flush=True)
        print("\n", file=control['h_log'],flush=True)
        print("\n", file=control['h_log'],flush=True)
        print("\n", file=control['h_log'],flush=True)
        control['iter_num_outer']=control['iter_num_outer']+1

    return None



def lqsgw_dmft(control,wan_hmat,imp):

    print("\n", file=control['h_log'],flush=True)
    print("\n", file=control['h_log'],flush=True)
    print("\n", file=control['h_log'],flush=True)

    print('*****   wannier  *****', file=control['h_log'],flush=True)
    if control['do_wannier']:

        check_wannier_function_input(control,wan_hmat)
        run_comwann(control, wan_hmat)
    comwann_postprocessing(control, wan_hmat)
    if control['do_wannier']:
        write_conv_wan(control)



    print('*****   Coulomb  *****', file=control['h_log'],flush=True)
    if control['do_coulomb']:    
        check_coulomb_input(control)
        run_comcoulomb(control,imp)
    comcoulomb_postprocessing(control,imp)
    if control['do_coulomb']:
        write_conv_coulomb(control,imp)


    print('*****   prepare dc  *****'    , file=control['h_log'],flush=True)


    if control['do_dc']:    
        prepare_initial_ef(control)
        generate_initial_transformation(control)        

        prepare_seed_dc_sig_and_wannier_dat(control,wan_hmat,imp)
        run_comlowh(control)

        prepare_dc(control,wan_hmat,imp) 
        run_dc(control,imp)

        cal_dc_diagonal(control)
        cal_zinv_m1_diagonal(control)

        generate_initial_self_energy(control,imp)
        write_conv_dc(control,imp)


    # control['iter_num_impurity']=control['iter_num_impurity']+1

    while (control['iter_num_impurity'] <= control['max_iter_num_impurity']):

        print('\n', file=control['h_log'],flush=True)
        print('*****   iter_num_impurity: ', str(control['iter_num_impurity']), '  *****', file=control['h_log'],flush=True)

        prepare_comlowh(control,wan_hmat,imp)
        run_comlowh(control)
        delta_causality=delta_postprocessing(control,imp)
        write_conv_delta(control,delta_causality)

        prepare_impurity_solver(control,wan_hmat,imp)
        run_impurity_solver(control,imp)


        control['iter_num_impurity']=control['iter_num_impurity']+1

    return None



def lqsgw_dmft_u_fixed(control,wan_hmat,imp):

    print("\n", file=control['h_log'],flush=True)
    print("\n", file=control['h_log'],flush=True)
    print("\n", file=control['h_log'],flush=True)

    print('*****   wannier  *****', file=control['h_log'],flush=True)
    if control['do_wannier']:

        check_wannier_function_input(control,wan_hmat)
        run_comwann(control, wan_hmat)
    comwann_postprocessing(control, wan_hmat)
    if control['do_wannier']:
        write_conv_wan(control)



    print('*****   Coulomb  *****', file=control['h_log'],flush=True)
    if control['do_coulomb']:    
        check_coulomb_input(control)
        run_comcoulomb(control,imp)
    comcoulomb_postprocessing(control,imp)
    if control['do_coulomb']:
        write_conv_coulomb(control,imp)


    print('*****   prepare dc  *****'    , file=control['h_log'],flush=True)


    if control['do_dc']:    
        prepare_initial_ef(control)
        generate_initial_transformation(control)        

        prepare_seed_dc_sig_and_wannier_dat(control,wan_hmat,imp)
        run_comlowh(control)

        prepare_dc(control,wan_hmat,imp) 
        run_dc(control,imp)

        cal_dc_diagonal_new(control)
        cal_zinv_m1_diagonal(control)

        generate_initial_self_energy(control,imp)
        write_conv_dc(control,imp)


    # control['iter_num_impurity']=control['iter_num_impurity']+1

    while (control['iter_num_impurity'] <= control['max_iter_num_impurity']):

        print('\n', file=control['h_log'],flush=True)
        print('*****   iter_num_impurity: ', str(control['iter_num_impurity']), '  *****', file=control['h_log'],flush=True)

        prepare_comlowh(control,wan_hmat,imp)
        run_comlowh(control)
        delta_causality=delta_postprocessing(control,imp)
        write_conv_delta(control,delta_causality)

        prepare_dc(control,wan_hmat,imp) 
        run_dc(control,imp)

        cal_dc_diagonal(control)
        cal_zinv_m1_diagonal(control)

        prepare_impurity_solver(control,wan_hmat,imp)
        run_impurity_solver(control,imp)

        control['iter_num_impurity']=control['iter_num_impurity']+1

    return None


if __name__ == '__main__':


    control,wan_hmat,imp=read_comdmft_ini()
    initial_file_directory_setup(control)    

    if (control['method'] == 'lda+dmft'):
        lda_dmft(control,wan_hmat,imp)
    elif (control['method'] == 'lqsgw+dmft'):
        lqsgw_dmft(control,wan_hmat,imp)
    # elif (control['method'] == 'lqsgw+dmft_u_fixed'):
    #     lqsgw_dmft_u_fixed(control,wan_hmat,imp)        
    close_h_log(control)



###### conv using tabulate

