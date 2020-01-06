import os, sys, subprocess
from pyglib.gutz import embedh
from pyglib.gsolver.gs_idmrg import create_idmrg_ctl_file


def driver_ridmrg(imp=1, lv2e=False):
    '''dmrg solver based on ITensor.
    '''
    if "-i" in sys.argv:
        imp = int(sys.argv[sys.argv.index("-i")+1])
    h1e, lambdac, daalpha, v2e = embedh.h5gen_embedh_spin_updn(imp=imp,
            lv2e=lv2e)
    embedh.wrt_text_rembed(h1e, lambdac, daalpha, v2e, shreshold=1.e-7)
    N = Npart = h1e.shape[0]
    create_idmrg_ctl_file(imp=imp, N=N, Npart=Npart)
    cmd = [os.environ["WIEN_GUTZ_ROOT2"]+"/exe_idmrg", str(imp)]
    with open("GDMRG_{}.LOG".format(imp), "w") as f:
        subprocess.call(cmd, stdout=f)
    dm, e_tot = embedh.get_res_idmrg(imp=imp)
    embedh.h5wrt_dm_sab_rc(dm, e_tot, imp=imp)


if __name__ == "__main__":
    driver_ridmrg(lv2e=True)
