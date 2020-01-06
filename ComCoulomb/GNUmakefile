#&&&&&&&&&&&&&&&&&&&&&&&&&&&&        UNIX MPICH MPI       &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

include ../arch.mk

# compfl = -debug -g -CB -check bounds -traceback -check uninit -fp-model precise

comp = $(PF90) $(FPPFLAGS) -c $(compfl)
link = $(PF90)  $(compfl)
LIB = $(LAPACK_LIB) 

PROG = ComCoulomb


.SUFFIXES: .mod .F .o



o_files = \
 etot_mod.o \
 atom_mod.o \
 manager_mod.o \
 parallel_mod.o \
 solid_mod.o \
 units_mod.o \
 vertex_mod.o \
 wannier90_mod.o \
 solid_scf_new.o \
 polarizability.o \
 screened_w.o \
 u_atomic.o \
 p_r_tau_solid.o \
 unport.o \
 etot_x_from_pw.o \
 pw_from_pb.o \
 tcgstb.o \
 boson_unpack_tau.o \
 check_w_sym.o \
 check_hermit.o \
 ws_k_nu_solid.o \
 output_u_f0f2f4f6.o \
 output_u_f0f2f4f6_rel.o \
 zone1_number.o \
 fft3.o \
 p_mm_r.o \
 sym_boson_rs.o \
 p_mm_q_from_r.o \
 p_im_r.o \
 p_ii_solid.o \
 sym_w_mt.o \
 sym_w_mi.o \
 gbs_number.o \
 p_im_q_from_r.o \
 ws_k_nu_solid_0.o \
 output_em1_gg0.o \
 eps_invers_gamma_prep2.o \
 output_w.o \
 wan_pb_rs.o \
 cft.o \
 lget.o \
 getlimj.o \
 g_rs_from_ks_mm.o \
 rotate.o \
 boson_pack_tau.o \
 g_rs_from_ks_mi.o \
 g_k_g_r1.o \
 g_k_g_r1_xi.o \
 g_rr_r_r1_star.o \
 s_overlap_m_b.o \
 sym_pro_w.o \
 sym_w_left.o \
 sym_a_0.o \
 sym_z_0.o \
 det3.o \
 sym_vmi_left.o \
 from_tau_to_nu_ba.o \
 pack_hermit.o \
 u_onsite.o \
 regularization.o \
 add_to_output_w.o \
 invers_z.o \
 add_to_output_e1.o \
 etot_c_from_pw.o \
 from_nu_to_tau_ba.o \
 output_boson_q_nu.o \
 g_x_tau.o \
 green0.o \
 ferm_unpack_tau.o \
 sym_val_left.o \
 prd1.o \
 prd2.o \
 from_tau_to_nu_cos_d.o \
 fint.o \
 sym_w.o \
 s_overlap.o \
 diagonalization_lapack.o \
 from_nu_to_tau_cos.o \
 sym_w_ii.o \
 solid.o \
 prod_bas.o \
 prod_b1.o \
 prod_b2.o \
 get_fi0_gw.o \
 get_tmm.o \
 get_fifi_rl.o \
 v_matr.o \
 v_im_r_from_k.o \
 size_shift_par.o \
 v_prn.o \
 dqdall.o \
 var_int.o \
 set_lattices.o \
 get_str_con.o \
 v_coul_mi.o  \
 v_coul_ii_gg.o \
 v_spectrum.o \
 v_mm_r_from_k.o \
 w_red_mi.o \
 v_coul_full.o \
 value_pb.o \
 sphharm.o \
 bessr.o \
 intalt.o \
 convx.o \
 lattc.o \
 rcnke0.o \
 cansl0.o \
 xxmste0.o \
 yymste.o \
 get_cfgt.o \
 diag_v_coul.o \
 n_overlap.o \
 sym_vmt0.o \
 fac2.o \
 lgen.o \
 dlmtor.o \
 erfxx.o \
 sym_pro.o \
 rotate11.o \
 fi0_full.o \
 z_wan_rr_mt.o \
 z_wan_rr_int.o \
 wannier_interface.o \
 read_gwl_bni.o \
 qft.o \
 beging.o \
 start_par_mb.o \
 start_par_solid.o \
 trans_pb_pw.o \
 correlated_atoms.o \
 read_ini.o \
 cal_weig_rc.o \
 read_ntau_bni.o \
 get_rr.o 

mod_files = \
 etot_mod.mod \
 atom_mod.mod \
 manager_mod.mod \
 parallel_mod.mod \
 solid_mod.mod \
 units_mod.mod \
 vertex_mod.mod \
 wannier90_mod.mod \
 lda_dmft_mod.mod

ifdef USE_HDF5
  o_files +=  read_gwl_bni_hdf5.o
  o_files +=  read_ntau_bni_hdf5.o
  o_files +=  hdf5io_mod.o
  o_files +=  read_wannier_dat_hdf5.o
  mod_files += hdf5io_mod.mod
endif

all: $(PROG)

mods: $(mod_files)

clean:
	rm -f *.o *.mod $(PROG)

.F.o:
	$(comp) $*.F

$(PROG) : $(mod_files) $(o_files)
	$(link) $(o_files) -o $(PROG) $(LIB)
	cp $(PROG) ../bin

.o.mod:
	echo $<

.c.o:
	mpicc -c $*.c

                                     
