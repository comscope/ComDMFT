import os, sys, subprocess
from pyglib.gutz import embedh


def driver_c2rspci_mott_onfly(imp=1):
    '''exe_spci_mott_onfly with complex embedding Hamiltonian to
    real transformation.
    '''
    if "-i" in sys.argv:
        imp = int(sys.argv[sys.argv.index("-i")+1])
    if os.path.isfile("EMBED_HAMIL_{}r.h5".format(imp)):
        lv2e = False
    else:
        lv2e = True
    # get embedding hamiltonian parameters.
    h1e, lambdac, daalpha, v2e = embedh.h5gen_embedh_spin_updn(imp=imp,
            lv2e=lv2e)
    # write real version of the embedding hamiltonian.
    embedh.h5wrt_rembed_hamil(h1e, lambdac, daalpha, v2e, imp=imp)

    cmd = [os.environ["WIEN_GUTZ_ROOT2"]+"/exe_rspci_mott_onfly", str(imp)]
    subprocess.call(cmd)

    # get results
    dm, e_mol = embedh.get_res_rspci_mott_onfly(imp=imp)
    embedh.h5wrt_dm_sab_rc(dm, e_mol, imp=imp)


if __name__ == "__main__":
    driver_c2rspci_mott_onfly()
