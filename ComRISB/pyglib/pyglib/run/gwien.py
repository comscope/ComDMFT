from __future__ import print_function
try:
    from builtins import range, zip
except:
    pass

import os, sys, glob, h5py, socket, shutil, time, re, \
        subprocess
import numpy as np
from collections import deque
from subprocess import Popen, PIPE
from pyglib.io.fio import file_exists
import pyglib.run.environ as env


def get_file_info(fname, unit, idmf, case, scratch, so, para, cmplx, _band,
        updn, dnup):
    '''help function to setup informations in def file.
    '''
    if 'in2' == fname:
        return [unit, "'{}.in2{}'".format(case, cmplx), "'old'",
                "'formatted'", 0]
    elif 'inso' == fname:
        return [unit, "'{}.inso'".format(case), "'unknown'", "'formatted'", 0]
    elif 'indmfl' == fname:
        return [unit, "'{}.indmfl'".format(case), "'old'", "'formatted'", 0]
    elif 'outputdmfupdn' == fname:
        return [unit, "'{}.outputdmf{}{}'".format(case, idmf, updn), \
                "'unknown'", "'formatted'", 0]
    elif 'in1c' == fname:
        return [unit, "'{}.in1c'".format(case), "'unknown'", "'formatted'", 0]
    elif 'vectorupdn' == fname:
        return [unit, "'{}/{}.vector{}{}{}'".format(scratch, case, so, \
                updn, para), "'unknown'","'unformatted'",9000]
    elif 'vectordnup' == fname:
        return [unit, "'{}/{}.vector{}{}{}'".format(scratch, case, so, \
                dnup, para), "'unknown'","'unformatted'",9000]
    elif 'klist' == fname:
        return [unit, "'{}.klist{}'".format(case, _band), "'old'", \
                "'formatted'", 0]
    elif 'kgen' == fname:
        return [unit, "'{}.kgen'".format(case), "'unknown'", "'formatted'", 0]
    elif 'vspupdn' == fname:
        return [unit, "'{}.vsp{}'".format(case, updn), "'old'", \
                "'formatted'", 0]
    elif 'vspdnup' == fname:
        return [unit, "'{}.vsp{}'".format(case, dnup), "'unknown'", \
                "'formatted'", 0]
    elif 'struct' == fname:
        return [unit, "'{}.struct'".format(case), "'old'", "'formatted'", 0]
    elif 'rotlm' == fname:
        return [unit, "'{}.rotlm'".format(case), "'unknown'", "'formatted'", 0]
    elif 'energysodum' == fname:
        if so == 'so':
            sodum = 'dum'
        else:
            sodum = dnup
        return [unit, "'{}.energy{}'".format(case, sodum), \
               "'unknown'", "'formatted'", 0]
    elif 'energyupdn' == fname:
        return [unit, "'{}.energy{}{}{}'".format(case, so, updn, para), \
                "'unknown'", "'formatted'", 0]
    elif 'energydnup' == fname:
        return [unit, "'{}.energy{}{}{}'".format(case, so, dnup, para), \
                "'unknown'", "'formatted'", 0]
    elif 'clmval' == fname:
        return [unit, "'{}.clmval{}'".format(case, updn), "'unknown'", \
                "'formatted'", 0]
    elif 'recprlist' == fname:
        return [unit, "'{}.recprlist'".format(case), "'unknown'", \
                "'formatted'", 9000]
    elif 'scf2updn' == fname:
        return [unit, "'{}.scf2{}'".format(case, updn), \
                "'unknown'", "'formatted'", 0]
    elif 'normupdn' == fname:
        if so == "so" and updn == "":
            _updn = "up"
        else:
            _updn = updn
        return [unit, "'{}.norm{}{}{}'".format(case, so, _updn, para), \
                "'unknown'", "'formatted'", 0]
    elif 'normdnup' == fname:
        return [unit, "'{}.norm{}{}{}'".format(case, so, dnup, para), \
                "'unknown'", "'formatted'", 0]
    else:
        raise ValueError('No matching file name {}!'.format(fname))


def fcreate_def_gwien(case, scratch='.', so='', para='', idmf='1', cmplx='',
        _band='', updn='', dnup='dn'):
    '''create gwien1/2.def file.
    '''
    fdef = open('gwien{}{}.def'.format(idmf,updn), 'w')
    if idmf == '1':
        fname_list = ['in2', 'inso', 'indmfl', 'outputdmfupdn', \
                'in1c', 'vectorupdn', 'vectordnup', 'klist', \
                'kgen', 'vspupdn', 'vspdnup', 'struct', \
                'rotlm', 'energydnup', 'energyupdn', 'normupdn', \
                'normdnup']
        unit_list = [3, 4, 5, 6, \
                7, 9, 10, 13, \
                14, 18, 19, 20, \
                22, 59, 60, 12, \
                11]
    elif idmf == '2':
        fname_list = ['in1c', 'inso', 'in2', 'outputdmfupdn', 'indmfl', \
                'clmval', 'vectorupdn', 'vectordnup', 'recprlist', 'kgen', \
                'vspupdn', 'struct', 'scf2updn', 'rotlm', 'energyupdn', \
                'normupdn', 'normdnup']
        unit_list = [3, 4, 5, 6, 7, \
                8, 9, 10, 13, 14, \
                18, 20, 21, 22, 30, \
                12, 11]

    for fname, unit in zip(fname_list, unit_list):
        fdef.write("{:3d}, {:<15s}, {:<10s}, {:<13s}, {:<4d}\n".format(\
                *get_file_info(fname, unit, idmf, case, scratch, so, \
                para, cmplx, _band, updn, dnup)))
    fdef.close()


def onestep(fday, case, exec_name, w_root, para="", so="", \
        band=None, updn=None):
    '''wien2k steps.
    '''
    time_start = time.strftime("%H:%M:%S")
    cmd = ['{}/x'.format(w_root), exec_name, '-f', case]
    if para != "":
        cmd.append(para)
    if band == '-band':
        cmd.append(band)
        if not os.path.isfile('EFLDA.INP'):
            shutil.copy2('EFLDA.OUT', 'EFLDA.INP')
    if updn in ["-up", "-dn"]:
        cmd.append(updn)
    if so == "so":
        cmd.extend(["-c", "-so"])

    print(' '.join(x for x in cmd))
    process = Popen(cmd, stdout=PIPE, stderr=PIPE)
    out, err = process.communicate()
    fday.write('>{:<10s} ({}) {}\n'.format(exec_name, time_start, out[:-1]))
    fday.flush()
    for f in glob.glob('{}.error*'.format(exec_name)):
        if os.path.getsize(f) > 0:
            print('error in {} from file: {}'.format(
                    f, open(f, 'r').readlines()))
            sys.exit(1)


def gonestep(fday, exec_name, mpi, updn=""):
    '''gwien1, CyGutz and gwien2 steps.
    '''
    time_start = time.strftime("%H:%M:%S")
    with open(':log', 'a') as f:
        f.write('{}>   {}\n'.format(time.strftime("%a %b %d %H:%M:%S %Z %Y"), \
                exec_name))

    cmd = ['/usr/bin/time']
    if mpi != '':
        cmd.extend(mpi)
    cmd.append('{}'.format(exec_name))
    if 'gwien' in exec_name:
        cmd.append('{}{}.def'.format(exec_name, updn))

    print(' '.join(x for x in cmd))
    process = Popen(cmd, stdout=PIPE, stderr=PIPE)
    out, err = process.communicate()
    with open('{}_info.out'.format(exec_name), 'w') as f:
        f.write(out)
    fday.write('>{:<10s} ({}) {}\n'.format(exec_name, time_start, \
            err.splitlines()[-2]))
    fday.flush()
    for f in glob.glob('{}.error*'.format(exec_name)):
        if os.path.getsize(f) > 0:
            print('error in {} from file: {}'.format(
                    f, open(f, 'r').readlines()))
            sys.exit(1)


def get_file_content(fname):
    if os.path.exists(fname):
        data = '\n------- {} --------\n'.format(fname)
        with open(fname, 'r') as f:
            data += f.read()
        return data
    else:
        return ''


def scf(case, spinpol):
    # scf file content
    if spinpol:
        f_list = ['{}.scf{}'.format(case, i) for i in ['0', \
                '1up', '1dn', 'so', '2up', '2dn', \
                '1s', '2s', 'cup', 'cdn']]
    else:
        f_list = ['{}.scf{}'.format(case, i) for i in ['0', \
                '1', 'so', '2', '1s', '2s', 'c']]

    data = ''.join(get_file_content(f) for f in f_list)

    with open('{}.scf'.format(case), 'a') as f:
        f.write(data)

    # files saved for mixing.
    if spinpol:
        f_list = ['clmsum', 'vspup', 'vspdn', 'vnsup', 'vnsdn', 'vrespsum',
                'clmdn', 'clmup']
    else:
        f_list = ['clmsum', 'vsp', 'vns', 'vrespsum']

    for i in f_list:
        name = '{}.{}'.format(case, i)
        if file_exists(name):
            shutil.copy2(name, '{}_old'.format(name))


def scfm(case):
    f_scf = '{}.scfm'.format(case)
    data = get_file_content(f_scf)

    with open('{}.scf'.format(case), 'a') as f:
        f.write(data)


def diff(fday, case, mix_dc, avg_dc, gskip):
    e_que = deque([], 2)
    with open('{}.scf'.format(case), 'r')  as f:
        for line in f:
            if ':DIS' in line:
                d_rho = float(line.split()[-1])
            if ':ENE' in line:
                e_que.append(float(line.split()[-1]))
    if len(e_que) == 2:
        d_etot = np.abs(e_que[1] - e_que[0])
    else:
        d_etot = 0.0

    dcv_err = 0.
    if not gskip:
        with h5py.File("GPARAM.h5", 'a') as f:
            ldc = f["/dc_mode"][0]
            if os.path.isfile("GDC_NELF_OUT.h5"):
                with h5py.File("GDC_NELF_OUT.h5", 'r') as fp:
                    nelf_list_inp = fp["/dc_nelf_list_inp"][()]
                    nelf_list_out = fp["/dc_nelf_list_out"][()]
                nelf_diff_list = nelf_list_out - nelf_list_inp
                nelf_list_mix = nelf_list_inp + mix_dc*nelf_diff_list
                if avg_dc:
                    valup = np.sum(nelf_list_mix[:,0])/nelf_list_mix.shape[0]
                    valdn = np.sum(nelf_list_mix[:,1])/nelf_list_mix.shape[0]
                    nelf_list_mix = [[valup,valdn] for x in nelf_list_inp]
                if ldc == 12:
                    if avg_dc:
                        dcv_err = np.sum(nelf_diff_list)/len(nelf_list_mix)
                    else:
                        dcv_err = np.max(np.abs(nelf_diff_list))
                    if '/dc_nelf_list' in f:
                        f["/dc_nelf_list"][()] = nelf_list_mix
                    else:
                        f["/dc_nelf_list"] = nelf_list_mix

    fday.write(':ENERGY convergence: {}\n'.format(d_etot))
    fday.write(':CHARGE convergence: {}\n'.format(d_rho))
    fday.write(':VDC convergence: {}\n'.format(dcv_err))
    return d_rho, d_etot, dcv_err


def processes_convert(so,updn):
    if not file_exists('.processes'):
        print('.processes file not present. It must be a serial run.')
        return
    lines = open('.processes').readlines()
    work = {}
    nkstart = 0
    for line in lines:
        data = line.split(':')
        if data[0].strip().isdigit():
            vecn = ["emmanuel" for i in range(6)]
            i, nkp, nprc = map(int,data[::2])
            if not so:
                fdef = open('{}lapw1_{}.def'.format(updn,i), 'r')
                for line in fdef:
                    data = line.split(',')
                    data0 = int(data[0])
                    if data0 == 10 or data0 == 11:
                        data0 = data0 % 10
                        m = re.search('.*[\'|\"](.*)_(\d+)', data[1])
                        assert m is not None, 'vector file to macth ' + \
                                ' lapw1.def not found!'
                        vecn[data0*2] = '{}_{}'.format(m.group(1), m.group(2))
                        vecn[data0*2+1] = '{}dn_{}'.format(m.group(1), \
                                 m.group(2))
                fdef.close()
            else:
                fdef = open('{}lapwso_{}.def'.format(updn,i), 'r')
                for line in fdef:
                    data = line.split(',')
                    if int(data[0])==42:
                        vecn[0]=data[1].split("'")[1]
                    elif int(data[0])==41:
                        vecn[1]=data[1].split("'")[1]
                    elif int(data[0])==52:
                        vecn[2]=data[1].split("'")[1]
                    elif int(data[0])==51:
                        vecn[3]=data[1].split("'")[1]
                    elif int(data[0])==46:
                        vecn[4]=data[1].split("'")[1]
                    elif int(data[0])==45:
                        vecn[5]=data[1].split("'")[1]
                fdef.close()

            if work.has_key(nprc):
                work[nprc].append((i, nkp, nkstart, vecn))
            else:
                work[nprc]=[(i, nkp, nkstart, vecn)]
            nkstart += nkp

    for prc in sorted(work.keys()):
        fo = open('_processes_{}'.format(prc-1), 'w')
        for (i, nkp, nkstart, vecn) in work[prc]:
            fo.write('{} {} {} "{}" "{}" "{}" "{}" "{}" "{}"\n'.format(\
                    i, nkp, nkstart, *vecn))


def create_gomp_file():
    '''
    Create GOMP.h5 file based on GMPI_X.h5 for openMP execution.
    '''
    with h5py.File('GMPI_0.h5', 'r') as f:
        num_procs = f["/nprocs"][0]
    nvec = 0
    kvec1 = []
    kvec2 = []
    for iproc in range(num_procs):
        with h5py.File('GMPI_' + str(iproc) + '.h5', 'r') as f:
            nvec += f["/nvec"][0]
            kvec = f["/KVEC"][()].T
            kvec1.append(kvec[0])
            if kvec.shape[1] == 2:
                kvec2.append(kvec[1])
    kvec = np.asarray(kvec1 + kvec2)

    with h5py.File('GOMP.h5', 'w') as f:
        f['/nvec'] = np.asarray([nvec])
        f['/KVEC'] = kvec.T


def run_gwien(nmaxiter=100, mix_dc=0.2, cc=1.e-3, ec=1.e-5, vc=1.e-2,
        startp='lapw0', endp='', band='', dos='',
        openmp=False, cygutz='CyGutz',
        pa_list=[], recycle_rl=True, avg_dc=True, spinpol=False,
        p_so=False, gskip=False):
    '''Driver for Wien2k + Gutzwiller-Slave-boson job.
    '''
    if '-s' in sys.argv:
        startp = sys.argv[sys.argv.index('-s') + 1]
    if '-e' in sys.argv:
        endp = sys.argv[sys.argv.index('-e') + 1]
    if '-cc' in sys.argv:
        cc = float(sys.argv[sys.argv.index('-cc') + 1])
    if '-ec' in sys.argv:
        ec = float(sys.argv[sys.argv.index('-ec') + 1])
    if '-vc' in sys.argv:
        vc = float(sys.argv[sys.argv.index('-vc') + 1])
    if '-n' in sys.argv:
        nmaxiter = int(sys.argv[sys.argv.index('-n') + 1])
    if '-omp' in sys.argv:
        openmp = True
        print('Using Open-MP instead of MPI of CyGutz.')
    if '-amix' in sys.argv:
        mix_dc = float(sys.argv[sys.argv.index('-amix') + 1])
    if '-band' in sys.argv:
        band ='-band'
    if '-dos' in sys.argv:
        dos = '-dos'
    if '-nrl' in sys.argv:
        recycle_rl = False
    if '-navg_dc' in sys.argv:
        avg_dc = False
    if "-sp" in sys.argv:
        spinpol = True
    if "-so" in sys.argv:
        p_so = True
    if "-dft" in sys.argv:
        gskip = True
    if band == '-band':
        _band = '_band'
        nmaxiter = 1
    else:
        _band = ''
    if band == '-band' or dos == '-dos':
        cygutz = 'CyGutzB'

    if len(sys.argv) > 1 and sys.argv[1] in ['-h', '--help']:
        help = '''
    The script is a wrapper to run Wien2k + Gutzwiller.
    It usually loops over the following steps:

        x lapw0   : computes LDA potential with current DFT+G-RISB charge
        x lapw1   : solves LDA eigenvalue equations
        [x lapwso]: second variational treatment of spin-orbit coupling
        x gwien1  : compute the local projector in the basis of DFT bands
        x cygutz  : solve the generic KS-Hubbard model using G-RISB
        x gwien2  : computes DFT+G-RISB valence charge
        x lcore   : computes DFT core charge
        x mixer   : mixes new charge density with the previous result

    The parameters with default values are as follows:

        name     default  inline-argument  help
        --------------------------------------------------------------------
        nmaxiter 100      -n 100           max charge mixing steps
        mix_dc   0.2      -amix            D.C. potential mxing param
        cc       1.e-3    -cc 1.e-3        charge density cutoff to exit
        ec       1.e-5    -ec 1.e-5        total energy cutoff to exit
        startp   'lapw0'  -s lapw0         start program
        endp     ''       -e ''            end program
        openmp   False    -omp             use openMP instead of openMPI
        rl       True     -nrl             start from previous GA solutions
        avg_dc   True     -navg_dc         average dc among atoms or not
        spinpol  False    -sp              spin-symmetry breaking at DFT level
        '''
        print(help)
        sys.exit(0)

    para = ''
    _para = ''
    if file_exists('.machines') :
        para = ' -p'
        _para = '_x'

    toclean = glob.glob('*.scf*') + glob.glob('*.error*') + \
            glob.glob('*.outputdmf?.*') + glob.glob('EMBED_HAMIL_RES*')
    for f in toclean:
        os.remove(f)

    struct_file = glob.glob('*.struct')
    if len(struct_file) != 1:
        raise ValueError('{} struct files present while only one must exist!'. \
                format(len(struct_file)))
    w_case = struct_file[0].split('.')[0]
    w_root = os.environ['WIENROOT']
    w_scratch = os.environ['SCRATCH']
    g_root = os.environ['WIEN_GUTZ_ROOT2']

    # infomation file
    fday = open(w_case + '.dayfile', 'w')
    fday.write('Calculating {} in {} \non {} with PID {}\n'.format(\
            w_case, os.getcwd(), socket.gethostname(), os.getpid()))

    print('calculation with spin-orbit = {}'.format(p_so))
    if p_so:
        so='so'
        cmplx = 'c'
    else:
        so = ''
        cmplx = ''

    # In addition, check in1c file
    if file_exists(w_case+'.in1c'):
        cmplx = 'c'

    f_mpi = 'mpi_prefix.dat'
    if os.path.isfile(f_mpi):
        with open(f_mpi, 'r') as f:
            mpi = f.readline().split()
        print('{} exists -- running in parallel mode.'.format(f_mpi))
        print(' '.join(x for x in mpi))
    else:
        if para != '':
            raise ValueError('missing mpi_prefix.dat with .machines present!')
        mpi = ''
        print('{} not available -- running in serial mode.'.format(f_mpi))

    if openmp:
        _mpi = ''
    else:
        _mpi = mpi

    if not gskip:
        p_list = ['gwien1', 'gwien2']
        p_list.append(cygutz)
        for pa in pa_list:
            if pa not in p_list:
                p_list.append(pa)
        for p in p_list:
            shutil.copy2(g_root+'/'+p, '.')

        # create gwien1/2.def files
        if spinpol:
            fcreate_def_gwien(w_case, scratch=w_scratch, \
                    so=so, para=_para, \
                    idmf='1', cmplx=cmplx, _band=_band, \
                    updn="up", dnup='dn')
            if not p_so:
                fcreate_def_gwien(w_case, scratch=w_scratch, \
                        so=so, para=_para,
                        idmf='1', cmplx=cmplx, _band=_band, \
                                updn="dn", dnup='up')

            fcreate_def_gwien(w_case, scratch=w_scratch, \
                    so=so, para=_para, \
                    idmf='2', cmplx=cmplx, _band=_band, \
                    updn="up", dnup='dn')
            fcreate_def_gwien(w_case, scratch=w_scratch, \
                    so=so, para=_para, \
                    idmf='2', cmplx=cmplx, _band=_band, \
                    updn="dn", dnup='up')
        else:
            fcreate_def_gwien(w_case, scratch=w_scratch, \
                    so=so, para=_para, \
                    idmf='1', cmplx=cmplx, _band=_band)
            fcreate_def_gwien(w_case, scratch=w_scratch, \
                    so=so, para=_para, \
                    idmf='2', cmplx=cmplx, _band=_band)

    # save SLURM_ related environment variables
    slurm_envs = env.get_env_dict(key="SLURM_")

    if nmaxiter > 0:
        if not os.path.isfile('{}.clmsum'.format(w_case)):
            err_msg = 'no {}.clmsum file found--necessary for lapw0!'.\
                    format(w_case)
            print(err_msg)
            fday.print(err_msg+'\n')
            sys.exit(1)
        for f in glob.glob('*.broyd*'):
            os.remove(f)

        fday.write('   start at {} with {} \n    1/{} to go.\n'.format(
                time.asctime(), startp, nmaxiter))

    # Major charge density loop
    for icycle in range(nmaxiter):

        # unset slurm environment variables for wien2k type run
        env.unset_environ(slurm_envs)

        if icycle > 0 or (startp in 'lapw0'):
            onestep(fday, w_case, 'lapw0', w_root, para=para)
        if icycle > 0 or (startp in 'lapw0 lapw1'):
            if spinpol:
                onestep(fday, w_case, 'lapw1', w_root, para=para, band=band, \
                        updn="-up")
                onestep(fday, w_case, 'lapw1', w_root, para=para, band=band, \
                        updn="-dn")
            else:
                onestep(fday, w_case, 'lapw1', w_root, para=para, band=band)
        if (icycle > 0 or (startp in 'lapw0 lapw1 lapwso')) and p_so:
            if spinpol:
                onestep(fday, w_case, 'lapwso', w_root, para=para, \
                        band=band, updn="-up")
            else:
                onestep(fday, w_case, 'lapwso', w_root, para=para, band=band)

        if icycle==0 and para != '':
            if spinpol:
                processes_convert(p_so, updn="up")
            else:
                processes_convert(p_so, updn="")

        #set slurm environment variables for mpi run
        env.set_environ(slurm_envs)

        if gskip:
            # run dft only
            if icycle > 0 or (startp in 'lapw0 lapw1 lapwso lapw2'):
                if spinpol:
                    onestep(fday, w_case, 'lapw2', w_root, para=para, \
                            updn="-up", so=so)
                    onestep(fday, w_case, 'lapw2', w_root, para=para, \
                            updn="-dn", so=so)
                else:
                    onestep(fday, w_case, 'lapw2', w_root, para=para, so=so)
        else:
            if icycle > 0 or (startp in 'lapw0 lapw1 lapwso gwien1'):
                if spinpol:
                    gonestep(fday, 'gwien1', mpi, updn="up")
                    if not p_so:
                        gonestep(fday, 'gwien1', mpi, updn="dn")
                else:
                    gonestep(fday, 'gwien1', mpi)
                if openmp:
                    create_gomp_file()
                elif os.path.isfile("GOMP.h5"):
                    os.remove("GOMP.h5")
            if endp == 'gwien1':
                sys.exit(0)

            if icycle > 0 or (startp in 'lapw0 lapw1 lapwso gwien1 CyGutz'):
                gonestep(fday, cygutz, _mpi)
            if band == '-band' or dos == '-dos':
                sys.exit(0)
            shutil.copy2('GUTZ.LOG', 'SAVE_GUTZ.LOG')
            if endp == 'CyGutz':
                sys.exit(0)
            if recycle_rl:
                shutil.copy2('WH_RL_OUT.h5', 'WH_RL_INP.h5')

            if spinpol:
                gonestep(fday, 'gwien2', mpi, updn="up")
                gonestep(fday, 'gwien2', mpi, updn="dn")
            else:
                gonestep(fday, 'gwien2', mpi)

            if endp == 'gwien2':
                sys.exit(0)


        # unset slurm environment variables for wien2k type run
        env.unset_environ(slurm_envs)

        if spinpol:
            onestep(fday, w_case, 'lcore', w_root, para='', updn="-up")
            onestep(fday, w_case, 'lcore', w_root, para='', updn="-dn")
        else:
            onestep(fday, w_case, 'lcore', w_root, para='')

        scf(w_case, spinpol)
        onestep(fday, w_case, 'mixer', w_root, para='')
        scfm(w_case)
        drho, dene, dvdc = diff(fday, w_case, mix_dc, avg_dc, gskip)

        if gskip:
            gerr = 0.
        else:
            with h5py.File('GLOG.h5', 'r') as f:
                gerr = f['/rl_maxerr'][0]

        print(('dc={:.1e}, cc={:.1e} -> {:.0e}, ec={:.1e} ' + \
                '-> {:.0e}, gc={:.1e} icycle={}').format(
                dvdc, drho, cc, dene, ec, gerr, icycle))
        if drho < cc and dene < ec and dvdc < vc:
            sys.exit(0)


def batch_init_ga(dir_template='./template'):
    '''Loop over all the directories to initialize CyGutz calculations
     -- actually, since the CyGutz input files remain the same for different
    volumes, it simply copy the input files in template directory to
    each folder.
    '''
    cwd = os.getcwd()+'/'
    for dname in glob.glob('V*'):
        os.chdir(dname+'/case')
        shutil.copy(cwd+'/'+dir_template+'/ginit.h5', './')
        shutil.copy(cwd+'/'+dir_template+'/GPARAM.h5', './')
        if os.path.isfile(cwd+'/'+dir_template+'/GESOLVER.h5'):
            shutil.copy(cwd+'/'+dir_template+'/GESOLVER.h5', './')
        shutil.copy(cwd+'/'+dir_template+'/case.indmfl', './')
        os.chdir(cwd)


def batch_init_mott(dir_template='./template'):
    '''Loop over all the directories to initialize CyGutz-Mott calculations
     -- actually, since the CyGutz input files remain the same for different
    volumes, it simply copy the input files in template directory to
    each folder.
    '''
    cwd = os.getcwd()+'/'
    for dname in glob.glob('V*'):
        os.chdir(dname+'/case')
        shutil.copy(cwd+'/'+dir_template+'/GMOTT.h5', './')
        os.chdir(cwd)


def batch_modify_ga_setup(args, nproc=1):
    '''Loop over all the directories to modify CyGutz set up file.
    '''
    cwd = os.getcwd()+'/'
    cmd = [os.environ['WIEN_GUTZ_ROOT2']+'/switch_gparam.py'] + args
    if '-p' in sys.argv:
        nproc = int(sys.argv[sys.argv.index('-p')+1])
    for i,dname in enumerate(glob.glob('V*')):
        os.chdir(dname+'/case')
        proc = subprocess.Popen(cmd)
        os.chdir(cwd)
        if (i+1) % nproc == 0:
            proc.communicate()


def batch_job_slurm(u, j, dir_template='./template', dir_work='./'):
    '''copy template/job.slurm file to each working directory and submit jobs.
    '''
    cwd = os.getcwd()+'/'
    jname = 'u{}j{}'.format(u, j)

    # command to modify u,j
    args = ['-unique_u_ev', u, '-unique_j_ev', j]
    cmd = [os.environ['WIEN_GUTZ_ROOT2']+'/switch_gparam.py'] + args

    if '-w' in sys.argv:
        dir_work = sys.argv[sys.argv.index('-w')+1]

    # command to submit job.
    cmd_s = ['qsub', './job.slurm']

    for dname in glob.glob('V*'):
        os.chdir(dname+'/case/'+dir_work)

        # get job.slurm file
        with open(cwd+'/'+dir_template+'/job.slurm', 'r') as fin:
            with open('./job.slurm', 'w') as fout:
                for line in fin:
                    fout.write(line.replace('VV', dname). \
                            replace('UJ', jname))
        # modify u,j
        proc = subprocess.Popen(cmd)
        proc.communicate()

        # submit
        proc = subprocess.Popen(cmd_s)
        os.chdir(cwd)


def run_ga(nproc=1):
    '''Loop over all the directories to run_ga using nproc processors.
    '''
    cmd = [os.environ['WIEN_GUTZ_ROOT2']+'/run_ga.py']
    cwd = os.getcwd()+'/'
    if '-p' in sys.argv:
        nproc = int(sys.argv[sys.argv.index('-p')+1])
    for i,dname in enumerate(glob.glob('V*')):
        os.chdir(dname+'/case')
        proc = subprocess.Popen(cmd)
        os.chdir(cwd)
        if (i+1) % nproc == 0:
            proc.communicate()


def batch_gsave(sdir='ldag', args=['-f']):
    '''Loop over all the directories to save_lapw.
    '''
    cmd = [os.environ['WIEN_GUTZ_ROOT2']+'/save_ldag', '-d'] + [sdir] + args
    cwd = os.getcwd()+'/'
    for dname in glob.glob('V*'):
        os.chdir(dname+'/case')
        subprocess.call(cmd)
        os.chdir(cwd)



if __name__=='__main__':
    fcreate_def_gwien('FeSb2', scratch='.', so='', para='',
            idmf='1', cmplx='', _band='', updn='', dnup='dn')
