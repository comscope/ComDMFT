import os, sys, subprocess
from pyglib.gutz import embedh


def create_idmrg_ctl_file(imp=1, N=14, Npart=14, Nsweeps=5, ecut=1.e-6,
        maxM=200):
    '''write idmrg control file.
    '''
    fname = "GDMRG_{}.CTL".format(imp)
    if os.path.isfile(fname):
        return
    with open(fname, "w") as f:
        f.write(
'''input
    {{
    N = {}
    Npart = {}
    ConserveSz = no
    quiet = yes
    Nsweeps = {}
    Ecut = {}
    maxM = {}
    }}
'''.format(N, Npart, Nsweeps, ecut, maxM))


def driver_idmrg(imp=1):
    '''dmrg solver based on ITensor.
    '''
    if "-i" in sys.argv:
        imp = int(sys.argv[sys.argv.index("-i")+1])
    if os.path.isfile("V2E_{}.INP".format(imp)):
        lv2e = False
    else:
        lv2e = True
    h1e, lambdac, daalpha, v2e = embedh.h5gen_embedh_spin_updn(
            imp=imp, lv2e=lv2e)
    embedh.wrt_text_cembed(h1e, lambdac, daalpha, v2e)
    N = Npart = h1e.shape[0]
    create_idmrg_ctl_file(imp=imp, N=N, Npart=Npart)
    cmd = [os.environ["WIEN_GUTZ_ROOT2"]+"/exe_idmrg", str(imp)]
    with open("GDMRG_{}.LOG".format(imp), "w") as f:
        subprocess.call(cmd, stdout=f)
    dm, e_tot = embedh.get_res_idmrg(imp=imp)
    embedh.h5wrt_dm_sab_cmplx(dm, e_tot, imp=imp)


if __name__ == "__main__":
    driver_idmrg()
