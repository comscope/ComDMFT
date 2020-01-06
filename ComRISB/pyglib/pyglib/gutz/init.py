from __future__ import print_function

import sys
from pyglib.gutz.init_magnet_conf import init_magnet_conf

def initialize():
    '''
    Initialization for the Gutzwiller-Slave-Boson solver.
    Store all the relevant information in GPARAM.5 file.
    '''
    log_file = open("init_ga.slog", 'w')
    print(" The executed commant line is:", file=log_file)
    print("    "+" ".join(sys.argv[:]), file=log_file)

    from pyglib.gutz.structure import get_gatoms, check_material
    from pyglib.gutz.gatom import h5calc_save_lrot
    # get the material
    material = get_gatoms()
    material.set_modify_mode()

    # Q & A
    from pyglib.gutz.usrqa import h5save_usr_qa_setup
    h5save_usr_qa_setup(material, log=log_file)

    # Setting the attributes to material
    material.h5set()

    # Set the self-energy structure
    material.set_SelfEnergy()

    # Set the Lie parameters for 3d-rotations
    material.set_LieParameters()

    # Set one-particle rotation matrix list
    material.set_one_particle_rotation_list()

    # Set Coulomb matrix list
    material.set_v2e_list()

    # Set S/L matrix vector in symmetry-adapted basis
    material.set_SL_vector_list()

    # check some info in log
    check_material(material, log_file)

    # save_gparam.h5
    from pyglib.gutz.ginput import save_gparam
    save_gparam(iso=material.iso, ispin=material.ispin,
            ityp_list=material.ityp_list,
            imap_list=material.imap_list, na2_list=material.na2_list,
            imix=material.gimix, iembeddiag=material.giembeddiag,
            sigma_list=material.sigma_list, v2e_list=material.v2e_list,
            sx_list=material.sx_list, sy_list=material.sy_list,
            sz_list=material.sz_list, lx_list=material.lx_list,
            ly_list=material.ly_list, lz_list=material.lz_list,
            utrans_list=material.utrans_list, ldc=material.ldc,
            u_avg_list=material.u_avg_list,
            j_avg_list=material.j_avg_list,
            nelf_list=material.nelf_list,
            rotations_list=material.rotations_list,
            lie_odd_params_list=material.Lie_Jodd_list,
            lie_even_params_list=material.Lie_Jeven_list,
            jgenerator_list=material.jgenerator_list,
            sp_rotations_list=material.sp_rotations_list,
            nval_bot_list=material.nval_bot_list,
            nval_top_list=material.nval_top_list)

    # local rotation matrix (l) in Hilbert space.
    h5calc_save_lrot(material)
    log_file.close()


if __name__ == '__main__':
    initialize()
