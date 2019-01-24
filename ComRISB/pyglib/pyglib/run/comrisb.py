from __future__ import print_function

import numpy as np
import os, sys, shutil, subprocess, glob, h5py, pickle, time
from tabulate import tabulate


def open_h_log(control):
    control['h_log'] = open('./cmd.log', 'w', 0)
    control['h_log'].write(
            "\n*********************************\n" +
            "             ComRISB\n" +
            "*********************************\n\n")

    control['h_conv'] = open('./convergence.log', 'w', 0)
    control['h_conv'].write(
            "\n*********************************\n" +
            "             ComRISB\n" +
            "*********************************\n\n")


def close_h_log(control):
    control['h_log'].close()
    control['h_conv'].close()


def read_comrisb_ini():
    vglobl = {}
    vlocal = {}
    execfile('comdmft.ini', vglobl, vlocal)
    control = vlocal['control']
    wan_hmat = vlocal['wan_hmat']
    imp = vlocal['imp']

    control['name'] = 'control'
    wan_hmat['name'] = 'wan_hmat'
    imp['name'] = 'imp'

    open_h_log(control)

    control['comsuitedir'] = os.environ.get('COMSUITE_BIN')
    control['conv_table'] = []

    # in control
    control['top_dir'] = os.path.abspath('./')
    control["spin_orbit"] = control.get("spin_orbit", False)
    check_key_in_string('impurity_problem', control)
    check_key_in_string('impurity_problem_equivalence', control)

    # g-risb specific
    control['crystal_field'] = control.get('crystal_field', True)
    control['full_orbital_polarization'] = \
            control.get('full_orbital_polarization', False)
    control['iembed_diag'] = control.get('iembed_diag', -3)
    control['lnewton'] = control.get('lnewton', 0)
    control['spin_polarization'] = control.get('spin_polarization', False)
    control['unit'] = control.get('unit', 'ev')

    control['max_iter_num_impurity'] = control.get('max_iter_num_impurity', 1)
    control['max_iter_num_outer'] = control.get('max_iter_num_outer', 50)
    control["initial_dft_dir"] = control.get('initial_dft_dir', "../dft/")
    control['initial_dft_dir'] = os.path.abspath(control['initial_dft_dir'])
    control['allfile'] = find_allfile(control['initial_dft_dir'])
    control['h_conv'].write('allfile = {}\n'.format(control['allfile']))
    control["mpi_prefix"] = control.get('mpi_prefix', "")
    control["mpi_prefix_wannier"] = control.get('mpi_prefix_wannier', \
            control["mpi_prefix"])
    control["mpi_prefix_lattice"] = control.get('mpi_prefix_lattice', \
            control["mpi_prefix"])

    if 'dc_mat_to_read' in control:
        control['dc_mat_to_read'] = os.path.abspath(control['dc_mat_to_read'])

    if 'additional_unitary_transform_impurity' not in control:
        control['additional_unitary_transform_impurity'] = 0

    if control['additional_unitary_transform_impurity'] != 0:
        check_key_in_string('metal_threshold', control)
        control['h_conv'].write(
                'impurity_orb_equivalence overwritten to be independent.\n')

    control['restart'] = control.get('restart', False)
    if "-c" in sys.argv:
        control['restart'] = True

    control['start_prog'] = control.get('start_prog', "initcopy")
    if "-s" in sys.argv:
        control['start_prog'] = sys.argv[sys.argv.index("-s")+1]
    control['end_prog'] = control.get('end_prog', "endconv")
    if "-e" in sys.argv:
        control['end_prog'] = sys.argv[sys.argv.index("-e")+1]

    if 'lattice_directory' not in control:
        control['lattice_directory'] = './lattice'
    if 'wannier_directory' not in control:
        control['wannier_directory'] = './wannier'
    control['lattice_directory'] = os.path.abspath(
            control['lattice_directory'])
    control['wannier_directory'] = os.path.abspath(
            control['wannier_directory'])

    if 'dc_directory' not in control:
        control['dc_directory'] = './dc'
    control['dc_directory'] = os.path.abspath(control['dc_directory'])
    if 'impurity_directory' not in control:
        control['impurity_directory'] = './impurity'
    control['impurity_directory'] = os.path.abspath(
            control['impurity_directory'])
    if 'lowh_directory' not in control:
        control['lowh_directory'] = './lowh'
    control['lowh_directory'] = os.path.abspath(control['lowh_directory'])
    if control['additional_unitary_transform_impurity'] != 0:
        check_key_in_string('initial_trans_basis', control)

    # in wan_hmat
    wan_hmat["kgrid"] = wan_hmat.get("kgrid", None)
    check_key_in_string('froz_win_min', wan_hmat)
    check_key_in_string('froz_win_max', wan_hmat)
    wan_hmat['dis_win_min'] = wan_hmat.get('dis_win_min',
            wan_hmat['froz_win_min'])
    wan_hmat['dis_win_max'] = wan_hmat.get('dis_win_max',
            wan_hmat['froz_win_max']+40.0)
    control['proj_win_min'] = control.get('proj_win_min',
            wan_hmat['froz_win_min'])
    control['proj_win_max'] = control.get('proj_win_max',
            wan_hmat['froz_win_max'])
    wan_hmat['num_iter'] = wan_hmat.get('num_iter', 0)
    wan_hmat['dis_num_iter'] = wan_hmat.get('dis_num_iter', 100)

    wan_hmat['cut_low']=wan_hmat.get('cut_low', 0.4)
    wan_hmat['cut_froz']=wan_hmat.get('cut_froz', 0.15)
    wan_hmat['cut_total']=wan_hmat.get('cut_total', 0.0)

    control['h_log'].write('top_dir {}\n'.format(control['top_dir']))
    control['h_log'].write('lattice_directory {}\n'.
            format(control['lattice_directory']))
    control['h_log'].write('wannier_directory {}\n'.
            format(control['wannier_directory']))
    control['h_log'].write('dc_directory {}\n'.format(\
            control['dc_directory']))
    control['h_log'].write('impurity_directory {}\n'.
            format(control['impurity_directory']))
    control['h_log'].write('lowh_directory {}\n'
            .format(control['lowh_directory']))

    return control, wan_hmat, imp


def find_impurity_wan(control, wan_hmat):
    num_wann = np.shape(wan_hmat['basis'])[0]
    control['impurity_wan'] = []
    for ip in range(np.shape(control['impurity_problem'])[0]):
        if control['spin_orbit']:
            if control['impurity_problem'][ip][1].lower() == 'f':
                control['impurity_wan'].append([0]*14)
                for iwan in range(num_wann):
                    if wan_hmat['basis'][iwan]['atom'] == \
                            control['impurity_problem'][ip][0] \
                            and wan_hmat['basis'][iwan]['l'] == 3:
                        if int(wan_hmat['basis'][iwan]['i']*2) == -1:
                            if int(wan_hmat['basis'][iwan]['m']*2) == -5:
                                control['impurity_wan'][ip][0] = \
                                        wan_hmat['basis'][iwan]['ind']
                            elif int(wan_hmat['basis'][iwan]['m']*2) == -3:
                                control['impurity_wan'][ip][1] = \
                                        wan_hmat['basis'][iwan]['ind']
                            elif int(wan_hmat['basis'][iwan]['m']*2) == -1:
                                control['impurity_wan'][ip][2] = \
                                        wan_hmat['basis'][iwan]['ind']
                            elif int(wan_hmat['basis'][iwan]['m']*2) == 1:
                                control['impurity_wan'][ip][3] = \
                                        wan_hmat['basis'][iwan]['ind']
                            elif int(wan_hmat['basis'][iwan]['m']*2) == 3:
                                control['impurity_wan'][ip][4] = \
                                        wan_hmat['basis'][iwan]['ind']
                            elif int(wan_hmat['basis'][iwan]['m']*2) == 5:
                                control['impurity_wan'][ip][5] = \
                                        wan_hmat['basis'][iwan]['ind']
                        elif int(wan_hmat['basis'][iwan]['i']*2) == 1:
                            if int(wan_hmat['basis'][iwan]['m']*2) == -7:
                                control['impurity_wan'][ip][6] = \
                                        wan_hmat['basis'][iwan]['ind']
                            elif int(wan_hmat['basis'][iwan]['m']*2) == -5:
                                control['impurity_wan'][ip][7] = \
                                        wan_hmat['basis'][iwan]['ind']
                            elif int(wan_hmat['basis'][iwan]['m']*2) == -3:
                                control['impurity_wan'][ip][8] = \
                                        wan_hmat['basis'][iwan]['ind']
                            elif int(wan_hmat['basis'][iwan]['m']*2) == -1:
                                control['impurity_wan'][ip][9] = \
                                        wan_hmat['basis'][iwan]['ind']
                            elif int(wan_hmat['basis'][iwan]['m']*2) == 1:
                                control['impurity_wan'][ip][10] = \
                                        wan_hmat['basis'][iwan]['ind']
                            elif int(wan_hmat['basis'][iwan]['m']*2) == 3:
                                control['impurity_wan'][ip][11] = \
                                        wan_hmat['basis'][iwan]['ind']
                            elif int(wan_hmat['basis'][iwan]['m']*2) == 5:
                                control['impurity_wan'][ip][12] = \
                                        wan_hmat['basis'][iwan]['ind']
                            elif int(wan_hmat['basis'][iwan]['m']*2) == 7:
                                control['impurity_wan'][ip][13] = \
                                        wan_hmat['basis'][iwan]['ind']
            if (control['impurity_wan'][ip].count(0) != 0):
                raise ValueError('Something wrong in find_impurity_wan.')
        else:
            if control['impurity_problem'][ip][1].lower() == 'd':
                control['impurity_wan'].append([0]*5)
                for iwan in range(num_wann):
                    if wan_hmat['basis'][iwan]['atom'] == \
                            control['impurity_problem'][ip][0]\
                            and wan_hmat['basis'][iwan]['l'] == 2:
                        if wan_hmat['basis'][iwan]['m'] == -2:
                            control['impurity_wan'][ip][0] = \
                                    wan_hmat['basis'][iwan]['ind']
                        elif wan_hmat['basis'][iwan]['m'] == -1:
                            control['impurity_wan'][ip][1] = \
                                    wan_hmat['basis'][iwan]['ind']
                        elif wan_hmat['basis'][iwan]['m'] == 0:
                            control['impurity_wan'][ip][2] = \
                                    wan_hmat['basis'][iwan]['ind']
                        elif wan_hmat['basis'][iwan]['m'] == 1:
                            control['impurity_wan'][ip][3] = \
                                    wan_hmat['basis'][iwan]['ind']
                        elif wan_hmat['basis'][iwan]['m'] == 2:
                            control['impurity_wan'][ip][4] = \
                                    wan_hmat['basis'][iwan]['ind']
            if (control['impurity_wan'][ip].count(0) != 0):
                raise ValueError(
                        'Something wrong in find_impurity_wan.')


def initial_file_directory_setup(control):
    tempdir = control['wannier_directory']
    if not os.path.isdir(tempdir):
        os.mkdir(tempdir)
    tempdir=control['lattice_directory']
    if not os.path.isdir(tempdir):
        os.mkdir(tempdir)
    tempdir=control['lowh_directory']
    if not os.path.isdir(tempdir):
        os.mkdir(tempdir)


def initial_lattice_directory_setup(control):
    os.chdir(control['lattice_directory'])
    files = glob.iglob(control['initial_dft_dir']+"/*.rst")
    for filename in files:
        shutil.copy(filename, './')
    files = glob.iglob(control['initial_dft_dir']+"/*el_density")
    for filename in files:
        shutil.copy(filename, './')
    if os.path.exists(control['initial_dft_dir']+'/kpath'):
        shutil.copy(control['initial_dft_dir']+'/kpath', './')
    if os.path.exists(control['initial_dft_dir']+'/ini'):
        shutil.copy(control['initial_dft_dir']+'/ini', './')
    else:
        raise ValueError("ini file not exists!")

    iter_string = '_'+str(control['iter_num_outer'])
    shutil.copy(control['initial_dft_dir']+'/'+control['allfile']+'.out',
            control['allfile']+iter_string+'.out')

    control['h_log'].write("initial dft directory setup done.\n")
    os.chdir(control['top_dir'])


def create_comwann_ini(control, wan_hmat):
    with open('comwann.ini', 'w') as f:
        f.write(control['lattice_directory']+'\n')
        f.write('dft\n')
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


def read_wan_hmat_basis(control):
    # in the wannier directory
    inip = np.loadtxt(control['wannier_directory']+'/wannier.inip')
    basis_info = []

    if (control['spin_orbit']):
        for ii in range(np.shape(inip)[0]):
            basis_info.append({'atom': int(inip[ii, 0]), \
                    'l': int(inip[ii, 1]), 'i': inip[ii, 2], \
                    'm': inip[ii, 3], 'xaxis': inip[ii, 4:7], \
                    'zaxis': inip[ii, 7:10], 'ind': ii})
    else:
        for ii in range(np.shape(inip)[0]):
            basis_info.append({'atom': int(inip[ii, 0]), \
                    'l': int(inip[ii, 1]), 'm': int(inip[ii, 2]), \
                    'xaxis': inip[ii, 3:6], 'zaxis': inip[ii, 6:9], \
                    'ind': ii})
    control['h_log'].write(\
            'reading wannier.inip to get basis information.\n')
    return basis_info


def check_key_in_string(key, dictionary):
    if key not in dictionary:
        raise ValueError('missing \''+key+'\' in '+dictionary['name'])


def labeling_file(filename, iter_string):
    dirname = os.path.abspath(os.path.dirname(filename))
    filenameonly = os.path.basename(filename)
    temp = filenameonly.split('.')
    shutil.copy(dirname+'/'+filenameonly, dirname+"/" +
                '.'.join(temp[0:-1])+iter_string+'.'+temp[-1])


def write_conv(control):
    os.chdir(control['lattice_directory'])
    iter_string = '_' + str(control['iter_num_outer'])
    with open('./dft'+iter_string+'.out', "r") as f:
        for line in f:
            if "charge density" in line:
                line = line.split()
                delta_rho=float(line[3])
            elif "ETOT" in line:
                line = line.split()
                etot = float(line[4])
        control['conv_table'].append([control['iter_num_outer'], \
                delta_rho, etot])
    with h5py.File(control['lowh_directory']+"/GLOG.h5", "r") as f:
        mu = f["/e_fermi"][0]
        err_risb = f["/rl_maxerr"][0]
        control['conv_table'][-1].extend([mu, err_risb])
        rmat2 = f["/BND_R"][()].swapaxes(1,2)
    # check quasiparticle weight
    zmin = 100.
    for rmat in rmat2:
        zmat = rmat.T.conj().dot(rmat)
        w = np.linalg.eigvalsh(zmat)
        zmin = min(zmin, np.amin(w))
    control['conv_table'][-1].append(zmin)
    with open(control['top_dir']+'/convergence.log', 'w') as outputfile:
        outputfile.write(tabulate(control['conv_table'], \
                headers=['i_outer', \
                'delta_rho', 'etot', "mu", 'err_risb', "min_z"], \
                numalign="right",  floatfmt=".8f"))
    os.chdir(control['top_dir'])


def check_wannier_function_input(control, wan_hmat):
    os.chdir(control['wannier_directory'])
    create_comwann_ini(control, wan_hmat)
    if os.path.exists(control['top_dir']+'/local_axis.dat'):
        shutil.copy(control['top_dir']+'/local_axis.dat', './')
    os.chdir(control['top_dir'])


def run_dft(control):
    os.chdir(control['lattice_directory'])
    iter_string = '_'+str(control['iter_num_outer'])
    cmd = control['mpi_prefix_lattice'] + ' ' + \
            control['comsuitedir'] + "/rspflapw.exe"
    with open(control['lattice_directory']+'/dft.out', 'w') as logfile:
        ret = subprocess.call(cmd, shell=True, stdout = logfile, \
                stderr = logfile)
        if ret != 0:
            raise ValueError("Error in dft. Check dft.out for error message.")
    allfile=control['allfile']
    labeling_file('./'+allfile+'.out',iter_string)
    shutil.move('./dft.out', './dft'+iter_string+'.out')
    control['h_log'].write("dft calculation done.\n")
    os.chdir(control['top_dir'])


def prepare_dft_input(control):
    shutil.copy(control['lowh_directory']+"/wannier_den_matrix.dat",
            control['lattice_directory'])
    control['h_log'].write("prepare_dft_input done.\n")


def hlog_time(f, prestr, endl=""):
    f.write("{} at {} {}".format(prestr, time.strftime("%d/%m/%Y %H:%M:%S"), \
            endl))


def wannier_run(control,wan_hmat,fullrun=True):
    os.chdir(control['wannier_directory'])
    if fullrun:
        cmd = control['mpi_prefix_wannier']+' '+\
                control['comsuitedir']+"/ComWann"
        control['h_log'].write(cmd+"\n")

        hlog_time(control['h_log'], "comwann start")
        with open(control['wannier_directory']+'/comwann.out', 'w') \
                as logfile:
            ret = subprocess.call(cmd, shell=True, stdout = logfile, \
                    stderr = logfile)
            if ret != 0:
                raise ValueError(
                        "Error in comwann. Check comwann.out or OUT.")
        hlog_time(control['h_log'], "end", endl="\n")

        iter_string = '_' + str(control['iter_num_outer'])
        shutil.move('./wannier.wout', './wannier'+iter_string+'.wout')

    wan_hmat['basis'] = read_wan_hmat_basis(control)
    find_impurity_wan(control, wan_hmat)
    control['h_log'].write("control['impurity_wan']: {}\n".format(\
            control['impurity_wan']))
    os.chdir(control['top_dir'])


def bool2yesno(boolean):
    if boolean:
        return "y"
    else:
        return "n"

def get_locrot_list(wan_hmat):
    iatm = -1
    lrot_list = []
    for basis in wan_hmat["basis"]:
        if basis["atom"] != iatm:
            x = np.asarray([float(e) for e in basis["xaxis"]])
            z = np.asarray([float(e) for e in basis["zaxis"]])
            y = np.cross(z, x)
            lrot_list.append(np.asarray([x,y,z]).T)
            iatm = basis["atom"]
    return lrot_list


def init_grisb(control, imp):
    with h5py.File("ginit.h5", "r") as f:
        symbols = f["/struct/symbols"][()]

    unique_df_list = []
    unique_f_list_ev = []
    unique_j_list_ev = []
    unique_u_list_ev = []
    unique_nf_list = []
    corr_atm_list = []
    corr_df_list = []
    for i,impurity in enumerate(control['impurity_problem']):
        corr_atm_list.append(impurity[0]-1) # to zero-base
        corr_df_list.append(impurity[1])
        if len(corr_atm_list) == 1 or \
                symbols[corr_atm_list[-1]] != symbols[corr_atm_list[-2]]:
            unique_df_list.append(impurity[1])
            f_list = np.zeros(4)
            for ifs,fs in enumerate(["f0", "f2", "f4", "f6"]):
                if fs in imp[str(i+1)]:
                    f_list[ifs] = imp[str(i+1)][fs]
                elif fs.upper() in imp[str(i+1)]:
                    f_list[ifs] = imp[str(i+1)][fs.upper()]
            unique_f_list_ev.append(f_list)
            unique_u_list_ev.append(f_list[0])
            if "nominal_n" in imp[str(i+1)]:
                nf = imp[str(i+1)]["nominal_n"]
                unique_nf_list.append([nf/2., nf/2.])
            if corr_df_list[-1].lower() == "d":
                j_hund = f_list[1]*(1.0 + 0.625)/14.0
            elif corr_df_list[-1].lower() == "f":
                j_hund = f_list[1]*(286.0 + 195.0*0.668 + \
                        250.0 *0.494)/6435.0
            elif corr_df_list[-1].lower() == "p":
                j_hund = f_list[1]/5.
            elif corr_df_list[-1].lower() == "s":
                j_hund = 0.
            else:
                raise ValueError("Not implemented for l = {}!".format(\
                        corr_df_list[-1]))
            unique_j_list_ev.append(j_hund)

        if len(corr_atm_list) > 1:
            # ensure atoms in ascendig order.
            assert corr_atm_list[-1] >= corr_atm_list[-2], \
                    "please specify impurities in asceding order!"

    if len(unique_nf_list) == 0:
        ldc = 12
    else:
        ldc = 2
    # Sanity check
    for i in range(1,len(corr_df_list)):
        if control["impurity_problem_equivalence"][i] == \
                control["impurity_problem_equivalence"][i-1]:
            assert corr_df_list[i] == corr_df_list[i-1], \
                    "correlated orbitals incompatible with equiv_list!"
    # setup unique_corr_symbol_list
    unique_corr_symbol_list = []
    for i in corr_atm_list:
        if symbols[i] not in unique_corr_symbol_list:
            unique_corr_symbol_list.append(symbols[i])
    idx_equivalent_atoms = []

    for i in range(len(symbols)):
        if i == 0:
            idx_equivalent_atoms.append(0)
            continue
        elif i in corr_atm_list:
            idx = corr_atm_list.index(i)
            if idx > 0:
                if control["impurity_problem_equivalence"][idx] \
                        == control["impurity_problem_equivalence"][idx-1]:
                    idx_equivalent_atoms.append(idx_equivalent_atoms[-1])
                    continue
        idx_equivalent_atoms.append(idx_equivalent_atoms[-1]+1)

    with h5py.File("ginit.h5", "a") as f:
        f["/usrqa/crystal_field"] = bool2yesno(\
                control["crystal_field"])
        f["/usrqa/full_orbital_polarization"] = bool2yesno(\
                control["full_orbital_polarization"])
        f["/usrqa/idx_equivalent_atoms"] = idx_equivalent_atoms
        f["/usrqa/iembeddiag"] = control["iembed_diag"]
        f["/usrqa/ldc"] = ldc
        f["/usrqa/lnewton"] = control['lnewton']
        f["/usrqa/spin_orbit_coup"] = bool2yesno(\
                control['spin_orbit'])
        f["/usrqa/spin_polarization"] = bool2yesno(\
                control['spin_polarization'])
        f["/usrqa/u_matrix_type"] = 3
        f["/usrqa/unique_corr_symbol_list"] = unique_corr_symbol_list
        f["/usrqa/unique_df_list"] = unique_df_list
        f["/usrqa/unique_j_list_ev"] = unique_j_list_ev
        f["/usrqa/unique_u_list_ev"] = unique_u_list_ev
        f["/usrqa/unique_f_list_ev"] = unique_f_list_ev
        if len(unique_nf_list) > 0:
            f["/usrqa/unique_nf_list"] = unique_nf_list
        f["/usrqa/unit"] =  control['unit']
    from pyglib.gutz.init import initialize as ginit
    ginit()


def gwannier_run(control, wan_hmat, imp, icycle):
    os.chdir(control['lowh_directory'])
    lrot_list = get_locrot_list(wan_hmat)
    params = {}
    params["corbs_list"] = control['impurity_wan']
    params["wpath"] = control['wannier_directory']
    params["lrot_list"] = lrot_list
    params["icycle"] = icycle
    if control['spin_orbit']:
        params["iso"] = 2
    params["lprefix"] = control['allfile']

    # prepare parameter file for mpirun
    with open("gwannier_params.pkl", "wb") as f:
        pickle.dump(params, f, protocol=pickle.HIGHEST_PROTOCOL)

    cmd = control['mpi_prefix_wannier'] + ' ' + \
            control['comsuitedir'] + "/gwannier.py"
    control['h_log'].write(cmd+"\n")
    hlog_time(control['h_log'], "gwannier start")
    with open(control['lowh_directory']+'/gwannier.out', 'w') as f:
        ret = subprocess.call(cmd, shell=True, stdout = f, \
                stderr = f)
    hlog_time(control['h_log'], "end", endl="\n")
    assert ret == 0, "Error in gwannier. Check gwannier.out."

    if icycle <= 1:
        init_grisb(control, imp)

    cmd = control['mpi_prefix_wannier'] + ' ' + \
            control['comsuitedir'] + "/CyGutz -r 0"
    control['h_log'].write(cmd+"\n")
    hlog_time(control['h_log'], "cygutz start")
    with open(control['lowh_directory']+'/grisb.out', 'w') as f:
        ret = subprocess.call(cmd, shell=True, stdout = f, \
                stderr = f)
    hlog_time(control['h_log'], "end", endl="\n")
    assert ret == 0, "Error in grisb. Check grisb.out."
    shutil.copy("./WH_RL_OUT.h5", "./WH_RL_INP.h5")
    shutil.copy("./GUTZ.LOG", "SAVE_GUTZ.LOG")

    cmd = control['mpi_prefix_wannier'] + ' ' + \
            control['comsuitedir'] + "/gwannden.py"
    control['h_log'].write(cmd+"\n")
    hlog_time(control['h_log'], "gwannden start")
    with open(control['lowh_directory']+'/gwannden.out', 'w') as f:
        ret = subprocess.call(cmd, shell=True, stdout = f, \
                stderr = f)
    hlog_time(control['h_log'], "end", endl="\n")
    assert ret == 0, "Error in gwannden. Check gwannden.out."
    os.chdir(control['top_dir'])


def find_allfile(dft_dir):
    files = glob.iglob(dft_dir+"/*.rst")
    for filename in files:
        temp = filename[:-4].split('/')[-1].split('_')
        if len(temp)==1 and temp[0]!='info':
            allfile = temp[0]
    try:
        return allfile
    except:
        raise ValueError("failed to locate *.rst files in {}!".format(\
                dft_dir))


def dft_risb(control, wan_hmat, imp):
    control['h_log'].write("\n\n")
    control['iter_num_outer'] = 1
    while control['iter_num_outer'] <= control['max_iter_num_outer']:
        control['h_log'].write(\
                "************************************************\n"+ \
                "iteration: {}\n".format(control['iter_num_outer'])+ \
                "************************************************\n")
        if control['iter_num_outer'] == 1 and \
                control['start_prog'] in ["initcopy"]:
            if not control['restart']:
                initial_lattice_directory_setup(control)

        fullrun = control['iter_num_outer'] > 1 or \
                control['start_prog'] in ["initcopy", "wannier"]
        check_wannier_function_input(control, wan_hmat)
        wannier_run(control, wan_hmat, fullrun)
        if control['iter_num_outer'] > 1 or \
                control['start_prog'] in ["initcopy", "wannier", "gwann"]:
            gwannier_run(control, wan_hmat, imp, control['iter_num_outer'])
        if control['iter_num_outer'] > 1 or \
                control['start_prog'] in ["initcopy", "wannier", \
                "gwann", "dft"]:
            prepare_dft_input(control)
            run_dft(control)
        if control['end_prog'] in ["dft"]:
            break
        write_conv(control)
        control['iter_num_outer']=control['iter_num_outer']+1



if __name__ == '__main__':
    control,wan_hmat,imp = read_comrisb_ini()
    initial_file_directory_setup(control)
    dft_risb(control, wan_hmat, imp)
    close_h_log(control)
