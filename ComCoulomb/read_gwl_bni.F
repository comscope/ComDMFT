      subroutine read_gwl_bni
      use manager_mod
      use atom_mod
      use solid_mod
      use parallel_mod
      use units_mod
      use vertex_mod
      use wannier90_mod
      use etot_mod
      implicit none
      include 'mpif.h'
      integer :: ierr
      integer :: i_len,k0,ind_k,ispin
      integer :: i_tau,ind_tau
      character*300 :: tau_ch
      integer :: i_nu,ind_nu
      character*300 :: nu_ch
      character*300 :: k0_ch
      character*1000 :: allfilet      
      integer :: sz_cg,sz_cgr,sz_indxcg,sz_jcg,sz_jcgr,
     $  sz_k_list1,sz_k_list2,sz_ksym_in_npnt
      logical :: is_cg,is_cgr,is_indxcg,is_jcg,is_jcgr,
     $  is_k_list,is_ksym_in_npnt
      logical :: is_augm,is_dr,is_gfun,is_gfund,is_h,is_ind_wf,
     $  is_indbasa,is_indfun0,is_io_lem,is_iob2,is_iopb,is_lf_isz,
     $  is_lfun,is_lfun_pb,is_lfunm,is_li_val,is_lim_pb_mt,is_lm_isz,
     $  is_lmb,is_lmpb,is_n_pbmt0,is_nrad,is_ntle,is_ntle_pb,is_r,
     $  is_r_atom,is_rl,is_smt,is_txtel,is_bett,is_et_sum_nu,
     $  is_ind_nu_int,is_indxcgr,is_tau_mesh,is_w_nu,is_w_nu_d,is_e_bnd,
     $  is_e_qp,is_g_sym_0,is_gbs,is_gbs_sym,is_i_kref,is_ifia,
     $  is_ig_pair,is_ind_k_star_back,is_ind_r_star,is_ind_r_star_back,
     $  is_indgb,is_indplw,is_indpw_gw,is_indpw_r,is_ip,is_ip_k,is_ip_r,
     $  is_iplf_gk,is_ir_ref,is_is,is_itshift_r,is_k_group,is_k_star,
     $  is_k_sym_0,is_list_sym_k,is_n_bnd,is_nbask,is_nhsort,is_nplw_r,
     $  is_nplwgw,is_nrr_red,is_num_r_star,is_num_sym_k,is_pnt,
     $  is_r0_pair,is_r_group,is_r_pnt,is_r_pnt_ind,is_shift,is_sovi,
     $  is_sovr,is_sym_label,is_tau,is_tshift,is_tshift_r,is_u,is_uj,
     $  is_wgt,is_igbs,is_iplf_bk
      logical :: is_ab_nu_from_tau,is_tau_from_nu,is_ev_bnd,
     $  is_q_qp,is_s_it_r0,is_v_ii,is_v_mm,is_w_rs_mi,is_z_bnd,is_ev_qp,
     $  is_z_qp
      if(maswrk) then
        open(1,file=adjustl(trim(rstpath))//'/gwl.bni',
     $    status='old',form='unformatted')
        read(1)limlb,limlb1,limlpb,max_e_pb,maxb,maxel,maxel_loc,
     $    maxlfpb,maxlfun,maxnrad,maxpb,maxpbr,maxntle,maxwf,max_rl,
     $    n_pbmt,n_pbmtm,n_pbtot,nfun,nlb2,nrmax,e_c_solid_pw,
     $    e_x_solid_pw,allfile,betta_t,chem_pot,chem_pot_qp,complex_ro,
     $    eps_pb,interp_nu_d,irel,itest,magn,maxylm,n_omega,n_nu,n_nu_d,
     $    n_nu_exa,n_tau,nrel,nrel_core,nspin,nspin_0,nspin_1,ubi,amega,
     $    gbas,k_star_max,max_r_star,maxwig,maxwig_j,metal,n_symb,natom,
     $    nbasmpw,nbndf,ndiv,ndtau,nkd,nkr,nplw_gw,nplw_rm,npnt,nqdiv,
     $    nr_full,nr_red,nrdiv,nrr_max,nrr_max_cor,nsort,numplw,par,
     $    r0b0,rb0,rbas,symmorph,v_gg_00,ngroup,maxplw,nplwdif,nbndtm,
     $    maxia,maxib,maxic,nplwbas,clight
        read(1)is_augm,is_dr,is_gfun,is_gfund,is_h,is_ind_wf,
     $    is_indbasa,is_indfun0,is_io_lem,is_iob2,is_iopb,is_lf_isz,
     $    is_lfun,is_lfun_pb,is_lfunm,is_li_val,is_lim_pb_mt,is_lm_isz,
     $    is_lmb,is_lmpb,is_n_pbmt0,is_nrad,is_ntle,is_ntle_pb,is_r,
     $    is_r_atom,is_rl,is_smt,is_txtel,is_bett,is_et_sum_nu,
     $    is_ind_nu_int,is_indxcgr,is_tau_mesh,is_w_nu,is_w_nu_d,
     $    is_e_bnd,is_e_qp,is_g_sym_0,is_gbs,is_gbs_sym,is_i_kref,
     $    is_ifia,is_ig_pair,is_ind_k_star_back,is_ind_r_star,
     $    is_ind_r_star_back,is_indgb,is_indplw,is_indpw_gw,is_indpw_r,
     $    is_ip,is_ip_k,is_ip_r,is_iplf_gk,is_ir_ref,is_is,is_itshift_r,
     $    is_k_group,is_k_star,is_k_sym_0,is_list_sym_k,
     $    is_n_bnd,is_nbask,is_nhsort,is_nplw_r,
     $    is_nplwgw,is_nrr_red,is_num_r_star,is_num_sym_k,is_pnt,
     $    is_r0_pair,is_r_group,is_r_pnt,is_r_pnt_ind,is_shift,is_sovi,
     $    is_sovr,is_sym_label,is_tau,is_tshift,is_tshift_r,is_u,is_uj,
     $    is_wgt,is_igbs,is_iplf_bk
      endif
      if(goparr) then
        call mpi_bcast(limlb,1,mpi_integer,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(limlb1,1,mpi_integer,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(limlpb,1,mpi_integer,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(max_e_pb,1,mpi_integer,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(maxb,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(maxel,1,mpi_integer,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(maxel_loc,1,mpi_integer,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(maxlfpb,1,mpi_integer,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(maxlfun,1,mpi_integer,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(maxnrad,1,mpi_integer,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(maxpb,1,mpi_integer,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(maxpbr,1,mpi_integer,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(maxntle,1,mpi_integer,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(maxwf,1,mpi_integer,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(max_rl,1,mpi_integer,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(n_pbmt,1,mpi_integer,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(n_pbmtm,1,mpi_integer,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(n_pbtot,1,mpi_integer,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(nfun,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(nlb2,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(nrmax,1,mpi_integer,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(e_c_solid_pw,1,mpi_double_precision,0,
     $    MPI_COMM_WORLD,mpierr)
        call mpi_bcast(e_x_solid_pw,1,mpi_double_precision,0,
     $    MPI_COMM_WORLD,mpierr)
        call mpi_bcast(allfile,72,mpi_character,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(betta_t,1,mpi_double_precision,0,
     $    MPI_COMM_WORLD,mpierr)
        call mpi_bcast(chem_pot,1,mpi_double_precision,0,
     $    MPI_COMM_WORLD,mpierr)
        call mpi_bcast(chem_pot_qp,1,mpi_double_precision,0,
     $    MPI_COMM_WORLD,mpierr)
        call mpi_bcast(complex_ro,1,mpi_logical,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(eps_pb,1,mpi_double_precision,0,
     $    MPI_COMM_WORLD,mpierr)
        call mpi_bcast(interp_nu_d,1,mpi_integer,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(irel,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(itest,1,mpi_integer,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(magn,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(maxylm,1,mpi_integer,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(n_omega,1,mpi_integer,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(n_nu,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(n_nu_d,1,mpi_integer,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(n_nu_exa,1,mpi_integer,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(n_tau,1,mpi_integer,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(nrel,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(nrel_core,1,mpi_integer,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(nspin,1,mpi_integer,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(nspin_0,1,mpi_integer,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(nspin_1,1,mpi_integer,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(ubi,3,mpi_character,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(amega,1,mpi_double_precision,0,
     $    MPI_COMM_WORLD,mpierr)
        call mpi_bcast(gbas,9,mpi_double_precision,0,
     $    MPI_COMM_WORLD,mpierr)
        call mpi_bcast(k_star_max,1,mpi_integer,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(max_r_star,1,mpi_integer,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(maxwig,1,mpi_integer,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(maxwig_j,1,mpi_integer,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(metal,1,mpi_logical,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(n_symb,1,mpi_integer,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(natom,1,mpi_integer,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(nbasmpw,1,mpi_integer,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(nbndf,1,mpi_integer,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(ndiv,3,mpi_integer,0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(ndtau,1,mpi_integer,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(nkd,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(nkr,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(nplw_gw,1,mpi_integer,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(nplw_rm,1,mpi_integer,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(npnt,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(nqdiv,1,mpi_integer,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(nr_full,1,mpi_integer,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(nr_red,1,mpi_integer,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(nrdiv,3,mpi_integer,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(nrr_max,1,mpi_integer,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(nrr_max_cor,1,mpi_integer,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(nsort,1,mpi_integer,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(numplw,1,mpi_integer,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(par,1,mpi_double_precision,0,
     $    MPI_COMM_WORLD,mpierr)
        call mpi_bcast(r0b0,9,mpi_double_precision,0,
     $    MPI_COMM_WORLD,mpierr)
        call mpi_bcast(rb0,9,mpi_double_precision,0,
     $    MPI_COMM_WORLD,mpierr)
        call mpi_bcast(rbas,9,mpi_double_precision,0,
     $    MPI_COMM_WORLD,mpierr)
        call mpi_bcast(symmorph,1,mpi_logical,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(v_gg_00,1,mpi_double_precision,0,
     $    MPI_COMM_WORLD,mpierr)
        call mpi_bcast(ngroup,1,mpi_integer,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(maxplw,1,mpi_integer,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(nplwdif,1,mpi_integer,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(nbndtm,1,mpi_integer,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(maxia,1,mpi_integer,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(maxib,1,mpi_integer,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(maxic,1,mpi_integer,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(nplwbas,1,mpi_integer,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(clight,1,mpi_double_precision,0,
     $    MPI_COMM_WORLD,mpierr)
      endif
      if(goparr) then
        call mpi_bcast(is_augm,1,mpi_logical,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(is_dr,1,mpi_logical,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(is_gfun,1,mpi_logical,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(is_gfund,1,mpi_logical,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(is_h,1,mpi_logical,0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(is_ind_wf,1,mpi_logical,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(is_indbasa,1,mpi_logical,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(is_indfun0,1,mpi_logical,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(is_io_lem,1,mpi_logical,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(is_iob2,1,mpi_logical,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(is_iopb,1,mpi_logical,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(is_lf_isz,1,mpi_logical,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(is_lfun,1,mpi_logical,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(is_lfun_pb,1,mpi_logical,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(is_lfunm,1,mpi_logical,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(is_li_val,1,mpi_logical,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(is_lim_pb_mt,1,mpi_logical,0,
     $    MPI_COMM_WORLD,mpierr)
        call mpi_bcast(is_lm_isz,1,mpi_logical,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(is_lmb,1,mpi_logical,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(is_lmpb,1,mpi_logical,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(is_n_pbmt0,1,mpi_logical,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(is_nrad,1,mpi_logical,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(is_ntle,1,mpi_logical,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(is_ntle_pb,1,mpi_logical,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(is_r,1,mpi_logical,0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(is_r_atom,1,mpi_logical,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(is_rl,1,mpi_logical,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(is_smt,1,mpi_logical,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(is_txtel,1,mpi_logical,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(is_bett,1,mpi_logical,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(is_et_sum_nu,1,mpi_logical,0,
     $    MPI_COMM_WORLD,mpierr)
        call mpi_bcast(is_ind_nu_int,1,mpi_logical,0,
     $    MPI_COMM_WORLD,mpierr)
        call mpi_bcast(is_indxcgr,1,mpi_logical,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(is_tau_mesh,1,mpi_logical,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(is_w_nu,1,mpi_logical,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(is_w_nu_d,1,mpi_logical,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(is_e_bnd,1,mpi_logical,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(is_e_qp,1,mpi_logical,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(is_g_sym_0,1,mpi_logical,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(is_gbs,1,mpi_logical,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(is_gbs_sym,1,mpi_logical,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(is_i_kref,1,mpi_logical,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(is_ifia,1,mpi_logical,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(is_ig_pair,1,mpi_logical,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(is_ind_k_star_back,1,mpi_logical,0,
     $    MPI_COMM_WORLD,mpierr)
        call mpi_bcast(is_ind_r_star,1,mpi_logical,0,
     $    MPI_COMM_WORLD,mpierr)
        call mpi_bcast(is_ind_r_star_back,1,mpi_logical,0,
     $    MPI_COMM_WORLD,mpierr)
        call mpi_bcast(is_indgb,1,mpi_logical,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(is_indplw,1,mpi_logical,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(is_indpw_gw,1,mpi_logical,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(is_indpw_r,1,mpi_logical,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(is_ip,1,mpi_logical,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(is_ip_k,1,mpi_logical,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(is_ip_r,1,mpi_logical,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(is_iplf_gk,1,mpi_logical,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(is_ir_ref,1,mpi_logical,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(is_is,1,mpi_logical,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(is_itshift_r,1,mpi_logical,0,
     $    MPI_COMM_WORLD,mpierr)
        call mpi_bcast(is_k_group,1,mpi_logical,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(is_k_star,1,mpi_logical,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(is_k_sym_0,1,mpi_logical,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(is_list_sym_k,1,mpi_logical,0,
     $    MPI_COMM_WORLD,mpierr)
        call mpi_bcast(is_n_bnd,1,mpi_logical,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(is_nbask,1,mpi_logical,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(is_nhsort,1,mpi_logical,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(is_nplw_r,1,mpi_logical,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(is_nplwgw,1,mpi_logical,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(is_nrr_red,1,mpi_logical,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(is_num_r_star,1,mpi_logical,0,
     $    MPI_COMM_WORLD,mpierr)
        call mpi_bcast(is_num_sym_k,1,mpi_logical,0,
     $    MPI_COMM_WORLD,mpierr)
        call mpi_bcast(is_pnt,1,mpi_logical,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(is_r0_pair,1,mpi_logical,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(is_r_group,1,mpi_logical,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(is_r_pnt,1,mpi_logical,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(is_r_pnt_ind,1,mpi_logical,0,
     $    MPI_COMM_WORLD,mpierr)
        call mpi_bcast(is_shift,1,mpi_logical,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(is_sovi,1,mpi_logical,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(is_sovr,1,mpi_logical,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(is_sym_label,1,mpi_logical,0,
     $    MPI_COMM_WORLD,mpierr)
        call mpi_bcast(is_tau,1,mpi_logical,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(is_tshift,1,mpi_logical,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(is_tshift_r,1,mpi_logical,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(is_u,1,mpi_logical,0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(is_uj,1,mpi_logical,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(is_wgt,1,mpi_logical,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(is_igbs,1,mpi_logical,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(is_iplf_bk,1,mpi_logical,0,MPI_COMM_WORLD,
     $    mpierr)
      endif
      if(is_augm .eq. .true.) then
        allocate(augm(maxntle,0:maxb,nsort))
      endif
      if(is_dr .eq. .true.) then
        allocate(dr(0:maxnrad,nsort))
      endif
      if(is_gfun .eq. .true.) then
        allocate(gfun(maxwf,nspin))
      endif
      if(is_gfund .eq. .true.) then
        allocate(gfund(maxwf,nspin))
      endif
      if(is_h .eq. .true.) then
        allocate(h(nsort))
      endif
      if(is_ind_wf .eq. .true.) then
        allocate(ind_wf(maxlfun,nsort))
      endif
      if(is_indbasa .eq. .true.) then
        allocate(indbasa(2,maxntle,limlb,nsort))
      endif
      if(is_indfun0 .eq. .true.) then
        allocate(indfun0(2,maxntle,nrel*maxb+1,nsort))
      endif
      if(is_io_lem .eq. .true.) then
        allocate(io_lem(natom))
      endif
      if(is_iob2 .eq. .true.) then
        allocate(iob2(natom))
      endif
      if(is_iopb .eq. .true.) then
        allocate(iopb(natom))
      endif
      if(is_lf_isz .eq. .true.) then
        allocate(lf_isz(maxel,nsort))
      endif
      if(is_lfun .eq. .true.) then
        allocate(lfun(nsort))
      endif
      if(is_lfun_pb .eq. .true.) then
        allocate(lfun_pb(nsort))
      endif
      if(is_lfunm .eq. .true.) then
        allocate(lfunm(nsort))
      endif
      if(is_li_val .eq. .true.) then
        allocate(li_val(maxlfun,nsort))
      endif
      if(is_lim_pb_mt .eq. .true.) then
        allocate(lim_pb_mt(0:maxpb,nsort))
      endif
      if(is_lm_isz .eq. .true.) then
        allocate(lm_isz(maxel,nsort))
      endif
      if(is_lmb .eq. .true.) then
        allocate(lmb(nsort))
      endif
      if(is_lmpb .eq. .true.) then
        allocate(lmpb(nsort))
      endif
      if(is_n_pbmt0 .eq. .true.) then
        allocate(n_pbmt0(nsort))
      endif
      if(is_nrad .eq. .true.) then
        allocate(nrad(nsort))
      endif
      if(is_ntle .eq. .true.) then
        allocate(ntle(0:maxb,nsort))
      endif
      if(is_ntle_pb .eq. .true.) then
        allocate(ntle_pb(0:maxpb,nsort))
      endif
      if(is_r .eq. .true.) then
        allocate(r(0:maxnrad,nsort))
      endif
      if(is_r_atom .eq. .true.) then
        allocate(r_atom(0:nrmax,nsort))
      endif
      if(is_rl .eq. .true.) then
        allocate(rl(0:maxnrad,0:max_rl,nsort))
      endif
      if(is_smt .eq. .true.) then
        allocate(smt(nsort))
      endif
      if(is_txtel .eq. .true.) then
        allocate(txtel(nsort))
      endif
      if(is_bett .eq. .true.) then
        allocate(bett((maxylm+1)**2))
      endif
      if(is_et_sum_nu .eq. .true.) then
        allocate(et_sum_nu(0:n_nu))
      endif
      if(is_ind_nu_int .eq. .true.) then
        allocate(ind_nu_int(n_nu))
      endif
      if(is_indxcgr .eq. .true.) then
        allocate(indxcgr(nrel_core*limlb1,nrel_core*limlb1,2,
     $    0:3*(magn-1)))
      endif
      if(is_tau_mesh .eq. .true.) then
        allocate(tau_mesh(0:n_tau))
      endif
      if(is_w_nu .eq. .true.) then
        allocate(w_nu(0:n_nu))
      endif
      if(is_w_nu_d .eq. .true.) then
        allocate(w_nu_d(0:n_nu_d))
      endif
      if(is_e_bnd .eq. .true.) then
        allocate(e_bnd(nbndf,npnt,nspin_0))
      endif
      if(is_e_qp .eq. .true.) then
        allocate(e_qp(nbndf,npnt,nspin_0))
      endif
      if(is_g_sym_0 .eq. .true.) then
        allocate(g_sym_0(nqdiv,ngroup))
      endif
      if(is_gbs .eq. .true.) then
        allocate(gbs(3,maxplw))
      endif
      if(is_gbs_sym .eq. .true.) then
        allocate(gbs_sym(nplwdif,ngroup))
      endif
      if(is_i_kref .eq. .true.) then
        allocate(i_kref(nqdiv))
      endif
      if(is_ifia .eq. .true.) then
        allocate(ifia(48))
      endif
      if(is_ig_pair .eq. .true.) then
        allocate(ig_pair(nqdiv,natom,natom))
      endif
      if(is_ind_k_star_back .eq. .true.) then
        allocate(ind_k_star_back(nqdiv))
      endif
      if(is_ind_r_star .eq. .true.) then
        allocate(ind_r_star(max_r_star,nr_red))
      endif
      if(is_ind_r_star_back .eq. .true.) then
        allocate(ind_r_star_back(nr_full))
      endif
      if(is_indgb .eq. .true.) then
        allocate(indgb(nbndtm,nqdiv))
      endif
      if(is_indplw .eq. .true.) then
        allocate(indplw(-maxia:maxia,-maxib:maxib,-maxic:maxic))
      endif
      if(is_indpw_gw .eq. .true.) then
        allocate(indpw_gw(nplw_gw,nqdiv))
      endif
      if(is_indpw_r .eq. .true.) then
        allocate(indpw_r(nplw_rm,nqdiv))
      endif
      if(is_ip .eq. .true.) then
        allocate(ip(natom,48))
      endif
      if(is_ip_k .eq. .true.) then
        allocate(ip_k(nqdiv,ngroup))
      endif
      if(is_ip_r .eq. .true.) then
        allocate(ip_r(nr_full,ngroup))
      endif
      if(is_iplf_gk .eq. .true.) then
        allocate(iplf_gk(nplwdif,nqdiv))
      endif
      if(is_ir_ref .eq. .true.) then
        allocate(ir_ref(nr_full))
      endif
      if(is_is .eq. .true.) then
        allocate(is(natom))
      endif
      if(is_itshift_r .eq. .true.) then
        allocate(itshift_r(3,nr_full,ngroup))
      endif
      if(is_k_group .eq. .true.) then
        allocate(k_group(nqdiv))
      endif
      if(is_k_star .eq. .true.) then
        allocate(k_star(npnt))
      endif
      if(is_k_sym_0 .eq. .true.) then
        allocate(k_sym_0(nqdiv,ngroup))
      endif
      if(is_list_sym_k .eq. .true.) then
        allocate(list_sym_k(ngroup,nqdiv))
      endif
      if(is_n_bnd .eq. .true.) then
        allocate(n_bnd(npnt,nspin_0),STAT=ierr)
      endif
      if(is_nbask .eq. .true.) then
        allocate(nbask(npnt),STAT=ierr)
      endif
      if(is_nhsort .eq. .true.) then
        allocate(nhsort(nsort))
      endif
      if(is_nplw_r .eq. .true.) then
        allocate(nplw_r(npnt))
      endif
      if(is_nplwgw .eq. .true.) then
        allocate(nplwgw(npnt))
      endif
      if(is_nrr_red .eq. .true.) then
        allocate(nrr_red(3,nqdiv*natom*natom))
      endif
      if(is_num_r_star .eq. .true.) then
        allocate(num_r_star(nr_red))
      endif
      if(is_num_sym_k .eq. .true.) then
        allocate(num_sym_k(nqdiv))
      endif
      if(is_pnt .eq. .true.) then
        allocate(pnt(3,nqdiv))
      endif
      if(is_r0_pair .eq. .true.) then
        allocate(r0_pair(nqdiv,natom,natom))
      endif
      if(is_r_group .eq. .true.) then
        allocate(r_group(nr_full))
      endif
      if(is_r_pnt .eq. .true.) then
        allocate(r_pnt(3,nr_full))
      endif
      if(is_r_pnt_ind .eq. .true.) then
        allocate(r_pnt_ind(nr_full))
      endif
      if(is_shift .eq. .true.) then
        allocate(shift(3,48))
      endif
      if(is_sovi .eq. .true.) then
        allocate(sovi(numplw))
      endif
      if(is_sovr .eq. .true.) then
        allocate(sovr(numplw))
      endif
      if(is_sym_label .eq. .true.) then
        allocate(sym_label(npnt))
      endif
      if(is_tau .eq. .true.) then
        allocate(tau(3,natom))
      endif
      if(is_tshift .eq. .true.) then
        allocate(tshift(3,natom,ngroup))
      endif
      if(is_tshift_r .eq. .true.) then
        allocate(tshift_r(3,nr_full,ngroup))
      endif
      if(is_u .eq. .true.) then
        allocate(u(maxwig,48))
      endif
      if(is_uj .eq. .true.) then
        allocate(uj(maxwig_j,48))
      endif
      if(is_wgt .eq. .true.) then
        allocate(wgt(nqdiv))
      endif
      if(is_igbs .eq. .true.) then
        allocate(igbs(3,maxplw))
      endif
      if(is_iplf_bk .eq. .true.) then
        allocate(iplf_bk(nplwbas,nqdiv))
      endif
      if(maswrk) then
        if(is_augm .eq. .true.) then
          read(1)augm
        endif
        if(is_dr .eq. .true.) then
          read(1)dr
        endif
        if(is_gfun .eq. .true.) then
          read(1)gfun
        endif
        if(is_gfund .eq. .true.) then
          read(1)gfund
        endif
        if(is_h .eq. .true.) then
          read(1)h
        endif
        if(is_ind_wf .eq. .true.) then
          read(1)ind_wf
        endif
        if(is_indbasa .eq. .true.) then
          read(1)indbasa
        endif
        if(is_indfun0 .eq. .true.) then
          read(1)indfun0
        endif
        if(is_io_lem .eq. .true.) then
          read(1)io_lem
        endif
        if(is_iob2 .eq. .true.) then
          read(1)iob2
        endif
        if(is_iopb .eq. .true.) then
          read(1)iopb
        endif
        if(is_lf_isz .eq. .true.) then
          read(1)lf_isz
        endif
        if(is_lfun .eq. .true.) then
          read(1)lfun
        endif
        if(is_lfun_pb .eq. .true.) then
          read(1)lfun_pb
        endif
        if(is_lfunm .eq. .true.) then
          read(1)lfunm
        endif
        if(is_li_val .eq. .true.) then
          read(1)li_val
        endif
        if(is_lim_pb_mt .eq. .true.) then
          read(1)lim_pb_mt
        endif
        if(is_lm_isz .eq. .true.) then
          read(1)lm_isz
        endif
        if(is_lmb .eq. .true.) then
          read(1)lmb
        endif
        if(is_lmpb .eq. .true.) then
          read(1)lmpb
        endif
        if(is_n_pbmt0 .eq. .true.) then
          read(1)n_pbmt0
        endif
        if(is_nrad .eq. .true.) then
          read(1)nrad
        endif
        if(is_ntle .eq. .true.) then
          read(1)ntle
        endif
        if(is_ntle_pb .eq. .true.) then
          read(1)ntle_pb
        endif
        if(is_r .eq. .true.) then
          read(1)r
        endif
        if(is_r_atom .eq. .true.) then
          read(1)r_atom
        endif
        if(is_rl .eq. .true.) then
          read(1)rl
        endif
        if(is_smt .eq. .true.) then
          read(1)smt
        endif
        if(is_txtel .eq. .true.) then
          read(1)txtel
        endif
        if(is_bett .eq. .true.) then
          read(1)bett
        endif
        if(is_et_sum_nu .eq. .true.) then
          read(1)et_sum_nu
        endif
        if(is_ind_nu_int .eq. .true.) then
          read(1)ind_nu_int
        endif
        if(is_indxcgr .eq. .true.) then
          read(1)indxcgr
        endif
        if(is_tau_mesh .eq. .true.) then
          read(1)tau_mesh
        endif
        if(is_w_nu .eq. .true.) then
          read(1)w_nu
        endif
        if(is_w_nu_d .eq. .true.) then
          read(1)w_nu_d
        endif
        if(is_e_bnd .eq. .true.) then
          read(1)e_bnd
        endif
        if(is_e_qp .eq. .true.) then
          read(1)e_qp
        endif
        if(is_g_sym_0 .eq. .true.) then
          read(1)g_sym_0
        endif
        if(is_gbs .eq. .true.) then
          read(1)gbs
        endif
        if(is_gbs_sym .eq. .true.) then
          read(1)gbs_sym
        endif
        if(is_i_kref .eq. .true.) then
          read(1)i_kref
        endif
        if(is_ifia .eq. .true.) then
          read(1)ifia
        endif
        if(is_ig_pair .eq. .true.) then
          read(1)ig_pair
        endif
        if(is_ind_k_star_back .eq. .true.) then
          read(1)ind_k_star_back
        endif
        if(is_ind_r_star .eq. .true.) then
          read(1)ind_r_star
        endif
        if(is_ind_r_star_back .eq. .true.) then
          read(1)ind_r_star_back
        endif
        if(is_indgb .eq. .true.) then
          read(1)indgb
        endif
        if(is_indplw .eq. .true.) then
          read(1)indplw
        endif
        if(is_indpw_gw .eq. .true.) then
          read(1)indpw_gw
        endif
        if(is_indpw_r .eq. .true.) then
          read(1)indpw_r
        endif
        if(is_ip .eq. .true.) then
          read(1)ip
        endif
        if(is_ip_k .eq. .true.) then
          read(1)ip_k
        endif
        if(is_ip_r .eq. .true.) then
          read(1)ip_r
        endif
        if(is_iplf_gk .eq. .true.) then
          read(1)iplf_gk
        endif
        if(is_ir_ref .eq. .true.) then
          read(1)ir_ref
        endif
        if(is_is .eq. .true.) then
          read(1)is
        endif
        if(is_itshift_r .eq. .true.) then
          read(1)itshift_r
        endif
        if(is_k_group .eq. .true.) then
          read(1)k_group
        endif
        if(is_k_star .eq. .true.) then
          read(1)k_star
        endif
        if(is_k_sym_0 .eq. .true.) then
          read(1)k_sym_0
        endif
        if(is_list_sym_k .eq. .true.) then
          read(1)list_sym_k
        endif
        if(is_n_bnd .eq. .true.) then
          read(1)n_bnd
        endif
        if(is_nbask .eq. .true.) then
          read(1)nbask
        endif
        if(is_nhsort .eq. .true.) then
          read(1)nhsort
        endif
        if(is_nplw_r .eq. .true.) then
          read(1)nplw_r
        endif
        if(is_nplwgw .eq. .true.) then
          read(1)nplwgw
        endif
        if(is_nrr_red .eq. .true.) then
          read(1)nrr_red
        endif
        if(is_num_r_star .eq. .true.) then
          read(1)num_r_star
        endif
        if(is_num_sym_k .eq. .true.) then
          read(1)num_sym_k
        endif
        if(is_pnt .eq. .true.) then
          read(1)pnt
        endif
        if(is_r0_pair .eq. .true.) then
          read(1)r0_pair
        endif
        if(is_r_group .eq. .true.) then
          read(1)r_group
        endif
        if(is_r_pnt .eq. .true.) then
          read(1)r_pnt
        endif
        if(is_r_pnt_ind .eq. .true.) then
          read(1)r_pnt_ind
        endif
        if(is_shift .eq. .true.) then
          read(1)shift
        endif
        if(is_sovi .eq. .true.) then
          read(1)sovi
        endif
        if(is_sovr .eq. .true.) then
          read(1)sovr
        endif
        if(is_sym_label .eq. .true.) then
          read(1)sym_label
        endif
        if(is_tau .eq. .true.) then
          read(1)tau
        endif
        if(is_tshift .eq. .true.) then
          read(1)tshift
        endif
        if(is_tshift_r .eq. .true.) then
          read(1)tshift_r
        endif
        if(is_u .eq. .true.) then
          read(1)u
        endif
        if(is_uj .eq. .true.) then
          read(1)uj
        endif
        if(is_wgt .eq. .true.) then
          read(1)wgt
        endif
        if(is_igbs .eq. .true.) then
          read(1)igbs
        endif
        if(is_iplf_bk .eq. .true.) then
          read(1)iplf_bk
        endif
      endif
      if(goparr) then
        if(is_augm .eq. .true.) then
          call mpi_bcast(augm,3*size(augm),mpi_character,0,
     $      MPI_COMM_WORLD,mpierr)
        endif
        if(is_dr .eq. .true.) then
          call mpi_bcast(dr,1*size(dr),mpi_double_precision,0,
     $      MPI_COMM_WORLD,mpierr)
        endif
        if(is_gfun .eq. .true.) then
          call mpi_bcast(gfun,1*size(gfun),mpi_double_precision,0,
     $      MPI_COMM_WORLD,mpierr)
        endif
        if(is_gfund .eq. .true.) then
          call mpi_bcast(gfund,1*size(gfund),mpi_double_precision,0,
     $      MPI_COMM_WORLD,mpierr)
        endif
        if(is_h .eq. .true.) then
          call mpi_bcast(h,1*size(h),mpi_double_precision,0,
     $      MPI_COMM_WORLD,mpierr)
        endif
        if(is_ind_wf .eq. .true.) then
          call mpi_bcast(ind_wf,1*size(ind_wf),mpi_integer,0,
     $      MPI_COMM_WORLD,mpierr)
        endif
        if(is_indbasa .eq. .true.) then
          call mpi_bcast(indbasa,1*size(indbasa),mpi_integer,0,
     $      MPI_COMM_WORLD,mpierr)
        endif
        if(is_indfun0 .eq. .true.) then
          call mpi_bcast(indfun0,1*size(indfun0),mpi_integer,0,
     $      MPI_COMM_WORLD,mpierr)
        endif
        if(is_io_lem .eq. .true.) then
          call mpi_bcast(io_lem,1*size(io_lem),mpi_integer,0,
     $      MPI_COMM_WORLD,mpierr)
        endif
        if(is_iob2 .eq. .true.) then
          call mpi_bcast(iob2,1*size(iob2),mpi_integer,0,
     $      MPI_COMM_WORLD,mpierr)
        endif
        if(is_iopb .eq. .true.) then
          call mpi_bcast(iopb,1*size(iopb),mpi_integer,0,
     $      MPI_COMM_WORLD,mpierr)
        endif
        if(is_lf_isz .eq. .true.) then
          call mpi_bcast(lf_isz,1*size(lf_isz),mpi_integer,0,
     $      MPI_COMM_WORLD,mpierr)
        endif
        if(is_lfun .eq. .true.) then
          call mpi_bcast(lfun,1*size(lfun),mpi_integer,0,
     $      MPI_COMM_WORLD,mpierr)
        endif
        if(is_lfun_pb .eq. .true.) then
          call mpi_bcast(lfun_pb,1*size(lfun_pb),mpi_integer,0,
     $      MPI_COMM_WORLD,mpierr)
        endif
        if(is_lfunm .eq. .true.) then
          call mpi_bcast(lfunm,1*size(lfunm),mpi_integer,0,
     $      MPI_COMM_WORLD,mpierr)
        endif
        if(is_li_val .eq. .true.) then
          call mpi_bcast(li_val,1*size(li_val),mpi_integer,0,
     $      MPI_COMM_WORLD,mpierr)
        endif
        if(is_lim_pb_mt .eq. .true.) then
          call mpi_bcast(lim_pb_mt,1*size(lim_pb_mt),mpi_integer,0,
     $      MPI_COMM_WORLD,mpierr)
        endif
        if(is_lm_isz .eq. .true.) then
          call mpi_bcast(lm_isz,1*size(lm_isz),mpi_integer,0,
     $      MPI_COMM_WORLD,mpierr)
        endif
        if(is_lmb .eq. .true.) then
          call mpi_bcast(lmb,1*size(lmb),mpi_integer,0,
     $      MPI_COMM_WORLD,mpierr)
        endif
        if(is_lmpb .eq. .true.) then
          call mpi_bcast(lmpb,1*size(lmpb),mpi_integer,0,
     $      MPI_COMM_WORLD,mpierr)
        endif
        if(is_n_pbmt0 .eq. .true.) then
          call mpi_bcast(n_pbmt0,1*size(n_pbmt0),mpi_integer,0,
     $      MPI_COMM_WORLD,mpierr)
        endif
        if(is_nrad .eq. .true.) then
          call mpi_bcast(nrad,1*size(nrad),mpi_integer,0,
     $      MPI_COMM_WORLD,mpierr)
        endif
        if(is_ntle .eq. .true.) then
          call mpi_bcast(ntle,1*size(ntle),mpi_integer,0,
     $      MPI_COMM_WORLD,mpierr)
        endif
        if(is_ntle_pb .eq. .true.) then
          call mpi_bcast(ntle_pb,1*size(ntle_pb),mpi_integer,0,
     $      MPI_COMM_WORLD,mpierr)
        endif
        if(is_r .eq. .true.) then
          call mpi_bcast(r,1*size(r),mpi_double_precision,0,
     $      MPI_COMM_WORLD,mpierr)
        endif
        if(is_r_atom .eq. .true.) then
          call mpi_bcast(r_atom,1*size(r_atom),mpi_double_precision,
     $      0,MPI_COMM_WORLD,mpierr)
        endif
        if(is_rl .eq. .true.) then
          call mpi_bcast(rl,1*size(rl),mpi_double_precision,0,
     $      MPI_COMM_WORLD,mpierr)
        endif
        if(is_smt .eq. .true.) then
          call mpi_bcast(smt,1*size(smt),mpi_double_precision,0,
     $      MPI_COMM_WORLD,mpierr)
        endif
        if(is_txtel .eq. .true.) then
          call mpi_bcast(txtel,4*size(txtel),mpi_character,0,
     $      MPI_COMM_WORLD,mpierr)
        endif
        if(is_bett .eq. .true.) then
          call mpi_bcast(bett,1*size(bett),mpi_double_precision,0,
     $      MPI_COMM_WORLD,mpierr)
        endif
        if(is_et_sum_nu .eq. .true.) then
          call mpi_bcast(et_sum_nu,1*size(et_sum_nu),
     $      mpi_double_precision,0,MPI_COMM_WORLD,mpierr)
        endif
        if(is_ind_nu_int .eq. .true.) then
          call mpi_bcast(ind_nu_int,1*size(ind_nu_int),mpi_integer,
     $      0,MPI_COMM_WORLD,mpierr)
        endif
        if(is_indxcgr .eq. .true.) then
          call mpi_bcast(indxcgr,1*size(indxcgr),mpi_integer,0,
     $      MPI_COMM_WORLD,mpierr)
        endif
        if(is_tau_mesh .eq. .true.) then
          call mpi_bcast(tau_mesh,1*size(tau_mesh),
     $      mpi_double_precision,0,MPI_COMM_WORLD,mpierr)
        endif
        if(is_w_nu .eq. .true.) then
          call mpi_bcast(w_nu,1*size(w_nu),mpi_double_precision,0,
     $      MPI_COMM_WORLD,mpierr)
        endif
        if(is_w_nu_d .eq. .true.) then
          call mpi_bcast(w_nu_d,1*size(w_nu_d),mpi_double_precision,
     $      0,MPI_COMM_WORLD,mpierr)
        endif
        if(is_e_bnd .eq. .true.) then
          call mpi_bcast(e_bnd,1*size(e_bnd),mpi_double_precision,0,
     $      MPI_COMM_WORLD,mpierr)
        endif
        if(is_e_qp .eq. .true.) then
          call mpi_bcast(e_qp,1*size(e_qp),mpi_double_precision,0,
     $      MPI_COMM_WORLD,mpierr)
        endif
        if(is_g_sym_0 .eq. .true.) then
          call mpi_bcast(g_sym_0,1*size(g_sym_0),mpi_integer,0,
     $      MPI_COMM_WORLD,mpierr)
        endif
        if(is_gbs .eq. .true.) then
          call mpi_bcast(gbs,1*size(gbs),mpi_double_precision,0,
     $      MPI_COMM_WORLD,mpierr)
        endif
        if(is_gbs_sym .eq. .true.) then
          call mpi_bcast(gbs_sym,1*size(gbs_sym),mpi_integer,0,
     $      MPI_COMM_WORLD,mpierr)
        endif
        if(is_i_kref .eq. .true.) then
          call mpi_bcast(i_kref,1*size(i_kref),mpi_integer,0,
     $      MPI_COMM_WORLD,mpierr)
        endif
        if(is_ifia .eq. .true.) then
          call mpi_bcast(ifia,1*size(ifia),mpi_integer,0,
     $      MPI_COMM_WORLD,mpierr)
        endif
        if(is_ig_pair .eq. .true.) then
          call mpi_bcast(ig_pair,1*size(ig_pair),mpi_integer,0,
     $      MPI_COMM_WORLD,mpierr)
        endif
        if(is_ind_k_star_back .eq. .true.) then
          call mpi_bcast(ind_k_star_back,1*size(ind_k_star_back),
     $      mpi_integer,0,MPI_COMM_WORLD,mpierr)
        endif
        if(is_ind_r_star .eq. .true.) then
          call mpi_bcast(ind_r_star,1*size(ind_r_star),mpi_integer,
     $      0,MPI_COMM_WORLD,mpierr)
        endif
        if(is_ind_r_star_back .eq. .true.) then
          call mpi_bcast(ind_r_star_back,1*size(ind_r_star_back),
     $      mpi_integer,0,MPI_COMM_WORLD,mpierr)
        endif
        if(is_indgb .eq. .true.) then
          call mpi_bcast(indgb,1*size(indgb),mpi_integer,0,
     $      MPI_COMM_WORLD,mpierr)
        endif
        if(is_indplw .eq. .true.) then
          call mpi_bcast(indplw,1*size(indplw),mpi_integer,0,
     $      MPI_COMM_WORLD,mpierr)
        endif
        if(is_indpw_gw .eq. .true.) then
          call mpi_bcast(indpw_gw,1*size(indpw_gw),mpi_integer,0,
     $      MPI_COMM_WORLD,mpierr)
        endif
        if(is_indpw_r .eq. .true.) then
          call mpi_bcast(indpw_r,1*size(indpw_r),mpi_integer,0,
     $      MPI_COMM_WORLD,mpierr)
        endif
        if(is_ip .eq. .true.) then
          call mpi_bcast(ip,1*size(ip),mpi_integer,0,MPI_COMM_WORLD,
     $      mpierr)
        endif
        if(is_ip_k .eq. .true.) then
          call mpi_bcast(ip_k,1*size(ip_k),mpi_integer,0,
     $      MPI_COMM_WORLD,mpierr)
        endif
        if(is_ip_r .eq. .true.) then
          call mpi_bcast(ip_r,1*size(ip_r),mpi_integer,0,
     $      MPI_COMM_WORLD,mpierr)
        endif
        if(is_iplf_gk .eq. .true.) then
          call mpi_bcast(iplf_gk,1*size(iplf_gk),mpi_integer,0,
     $      MPI_COMM_WORLD,mpierr)
        endif
        if(is_ir_ref .eq. .true.) then
          call mpi_bcast(ir_ref,1*size(ir_ref),mpi_integer,0,
     $      MPI_COMM_WORLD,mpierr)
        endif
        if(is_is .eq. .true.) then
          call mpi_bcast(is,1*size(is),mpi_integer,0,MPI_COMM_WORLD,
     $      mpierr)
        endif
        if(is_itshift_r .eq. .true.) then
          call mpi_bcast(itshift_r,1*size(itshift_r),mpi_integer,0,
     $      MPI_COMM_WORLD,mpierr)
        endif
        if(is_k_group .eq. .true.) then
          call mpi_bcast(k_group,1*size(k_group),mpi_integer,0,
     $      MPI_COMM_WORLD,mpierr)
        endif
        if(is_k_star .eq. .true.) then
          call mpi_bcast(k_star,1*size(k_star),mpi_integer,0,
     $      MPI_COMM_WORLD,mpierr)
        endif
        if(is_k_sym_0 .eq. .true.) then
          call mpi_bcast(k_sym_0,1*size(k_sym_0),mpi_integer,0,
     $      MPI_COMM_WORLD,mpierr)
        endif
        if(is_list_sym_k .eq. .true.) then
          call mpi_bcast(list_sym_k,1*size(list_sym_k),mpi_integer,
     $      0,MPI_COMM_WORLD,mpierr)
        endif
        if(is_n_bnd .eq. .true.) then
          call mpi_bcast(n_bnd,1*size(n_bnd),mpi_integer,0,
     $      MPI_COMM_WORLD,mpierr)
        endif
        if(is_nbask .eq. .true.) then
          call mpi_bcast(nbask,1*size(nbask),mpi_integer,0,
     $      MPI_COMM_WORLD,mpierr)
        endif
        if(is_nhsort .eq. .true.) then
          call mpi_bcast(nhsort,1*size(nhsort),mpi_integer,0,
     $      MPI_COMM_WORLD,mpierr)
        endif
        if(is_nplw_r .eq. .true.) then
          call mpi_bcast(nplw_r,1*size(nplw_r),mpi_integer,0,
     $      MPI_COMM_WORLD,mpierr)
        endif
        if(is_nplwgw .eq. .true.) then
          call mpi_bcast(nplwgw,1*size(nplwgw),mpi_integer,0,
     $      MPI_COMM_WORLD,mpierr)
        endif
        if(is_nrr_red .eq. .true.) then
          call mpi_bcast(nrr_red,1*size(nrr_red),mpi_integer,0,
     $      MPI_COMM_WORLD,mpierr)
        endif
        if(is_num_r_star .eq. .true.) then
          call mpi_bcast(num_r_star,1*size(num_r_star),mpi_integer,
     $      0,MPI_COMM_WORLD,mpierr)
        endif
        if(is_num_sym_k .eq. .true.) then
          call mpi_bcast(num_sym_k,1*size(num_sym_k),mpi_integer,0,
     $      MPI_COMM_WORLD,mpierr)
        endif
        if(is_pnt .eq. .true.) then
          call mpi_bcast(pnt,1*size(pnt),mpi_double_precision,0,
     $      MPI_COMM_WORLD,mpierr)
        endif
        if(is_r0_pair .eq. .true.) then
          call mpi_bcast(r0_pair,1*size(r0_pair),mpi_integer,0,
     $      MPI_COMM_WORLD,mpierr)
        endif
        if(is_r_group .eq. .true.) then
          call mpi_bcast(r_group,1*size(r_group),mpi_integer,0,
     $      MPI_COMM_WORLD,mpierr)
        endif
        if(is_r_pnt .eq. .true.) then
          call mpi_bcast(r_pnt,1*size(r_pnt),mpi_integer,0,
     $      MPI_COMM_WORLD,mpierr)
        endif
        if(is_r_pnt_ind .eq. .true.) then
          call mpi_bcast(r_pnt_ind,1*size(r_pnt_ind),mpi_integer,0,
     $      MPI_COMM_WORLD,mpierr)
        endif
        if(is_shift .eq. .true.) then
          call mpi_bcast(shift,1*size(shift),mpi_double_precision,0,
     $      MPI_COMM_WORLD,mpierr)
        endif
        if(is_sovi .eq. .true.) then
          call mpi_bcast(sovi,1*size(sovi),mpi_double_precision,0,
     $      MPI_COMM_WORLD,mpierr)
        endif
        if(is_sovr .eq. .true.) then
          call mpi_bcast(sovr,1*size(sovr),mpi_double_precision,0,
     $      MPI_COMM_WORLD,mpierr)
        endif
        if(is_sym_label .eq. .true.) then
          call mpi_bcast(sym_label,1*size(sym_label),mpi_character,
     $      0,MPI_COMM_WORLD,mpierr)
        endif
        if(is_tau .eq. .true.) then
          call mpi_bcast(tau,1*size(tau),mpi_double_precision,0,
     $      MPI_COMM_WORLD,mpierr)
        endif
        if(is_tshift .eq. .true.) then
          call mpi_bcast(tshift,1*size(tshift),mpi_double_precision,
     $      0,MPI_COMM_WORLD,mpierr)
        endif
        if(is_tshift_r .eq. .true.) then
          call mpi_bcast(tshift_r,1*size(tshift_r),
     $      mpi_double_precision,0,MPI_COMM_WORLD,mpierr)
        endif
        if(is_u .eq. .true.) then
          call mpi_bcast(u,1*size(u),mpi_double_precision,0,
     $      MPI_COMM_WORLD,mpierr)
        endif
        if(is_uj .eq. .true.) then
          call mpi_bcast(uj,1*size(uj),mpi_double_complex,0,
     $      MPI_COMM_WORLD,mpierr)
        endif
        if(is_wgt .eq. .true.) then
          call mpi_bcast(wgt,1*size(wgt),mpi_double_precision,0,
     $      MPI_COMM_WORLD,mpierr)
        endif
        if(is_igbs .eq. .true.) then
          call mpi_bcast(igbs,1*size(igbs),mpi_integer,0,
     $      MPI_COMM_WORLD,mpierr)
        endif
        if(is_iplf_bk .eq. .true.) then
          call mpi_bcast(iplf_bk,1*size(iplf_bk),mpi_integer,0,
     $      MPI_COMM_WORLD,mpierr)
        endif
      endif
      if(maswrk) then
        read(1)sz_cg,sz_cgr,sz_indxcg,sz_jcg,sz_jcgr,sz_k_list1,
     $    sz_k_list2,sz_ksym_in_npnt
        read(1)is_cg,is_cgr,is_indxcg,is_jcg,is_jcgr,is_k_list,
     $    is_ksym_in_npnt
      endif
      if(goparr) then
        call mpi_bcast(sz_cg,1,mpi_integer,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(sz_cgr,1,mpi_integer,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(sz_indxcg,1,mpi_integer,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(sz_jcg,1,mpi_integer,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(sz_jcgr,1,mpi_integer,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(sz_k_list1,1,mpi_integer,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(sz_k_list2,1,mpi_integer,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(sz_ksym_in_npnt,1,mpi_integer,0,
     $    MPI_COMM_WORLD,mpierr)
      endif
      if(goparr) then
        call mpi_bcast(is_cg,1,mpi_logical,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(is_cgr,1,mpi_logical,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(is_indxcg,1,mpi_logical,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(is_jcg,1,mpi_logical,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(is_jcgr,1,mpi_logical,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(is_k_list,1,mpi_logical,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(is_ksym_in_npnt,1,mpi_logical,0,
     $    MPI_COMM_WORLD,mpierr)
      endif
      if(is_cg .eq. .true.) then
        allocate(cg(sz_cg))
      endif
      if(is_cgr .eq. .true.) then
        allocate(cgr(sz_cgr))
      endif
      if(is_indxcg .eq. .true.) then
        allocate(indxcg(sz_indxcg))
      endif
      if(is_jcg .eq. .true.) then
        allocate(jcg(sz_jcg))
      endif
      if(is_jcgr .eq. .true.) then
        allocate(jcgr(sz_jcgr))
      endif
      if(is_k_list .eq. .true.) then
        allocate(k_list(sz_k_list1,sz_k_list2))
      endif
      if(is_ksym_in_npnt .eq. .true.) then
        allocate(ksym_in_npnt(sz_ksym_in_npnt))
      endif
      if(maswrk) then
        if(is_cg .eq. .true.) then
          read(1)cg
        endif
        if(is_cgr .eq. .true.) then
          read(1)cgr
        endif
        if(is_indxcg .eq. .true.) then
          read(1)indxcg
        endif
        if(is_jcg .eq. .true.) then
          read(1)jcg
        endif
        if(is_jcgr .eq. .true.) then
          read(1)jcgr
        endif
        if(is_k_list .eq. .true.) then
          read(1)k_list
        endif
        if(is_ksym_in_npnt .eq. .true.) then
          read(1)ksym_in_npnt
        endif
      endif
      if(goparr) then
        if(is_cg .eq. .true.) then
          call mpi_bcast(cg,1*size(cg),mpi_double_precision,0,
     $      MPI_COMM_WORLD,mpierr)
        endif
        if(is_cgr .eq. .true.) then
          call mpi_bcast(cgr,1*size(cgr),mpi_double_complex,0,
     $      MPI_COMM_WORLD,mpierr)
        endif
        if(is_indxcg .eq. .true.) then
          call mpi_bcast(indxcg,1*size(indxcg),mpi_integer,0,
     $      MPI_COMM_WORLD,mpierr)
        endif
        if(is_jcg .eq. .true.) then
          call mpi_bcast(jcg,1*size(jcg),mpi_integer,0,
     $      MPI_COMM_WORLD,mpierr)
        endif
        if(is_jcgr .eq. .true.) then
          call mpi_bcast(jcgr,1*size(jcgr),mpi_integer,0,
     $      MPI_COMM_WORLD,mpierr)
        endif
        if(is_k_list .eq. .true.) then
          call mpi_bcast(k_list,1*size(k_list),mpi_integer,0,
     $      MPI_COMM_WORLD,mpierr)
        endif
        if(is_ksym_in_npnt .eq. .true.) then
          call mpi_bcast(ksym_in_npnt,1*size(ksym_in_npnt),
     $      mpi_integer,0,MPI_COMM_WORLD,mpierr)
        endif
      endif
      call start_par_solid
      if(maswrk) then
        read(1)is_ab_nu_from_tau,is_tau_from_nu,is_ev_bnd,is_q_qp,
     $    is_s_it_r0,is_v_ii,is_v_mm,is_w_rs_mi,is_z_bnd,
     $    is_ev_qp,is_z_qp
      endif
      if(goparr) then
        call mpi_bcast(is_ab_nu_from_tau,1,mpi_logical,0,
     $    MPI_COMM_WORLD,mpierr)
        call mpi_bcast(is_tau_from_nu,1,mpi_logical,0,
     $    MPI_COMM_WORLD,mpierr)
        call mpi_bcast(is_ev_bnd,1,mpi_logical,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(is_q_qp,1,mpi_logical,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(is_s_it_r0,1,mpi_logical,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(is_v_ii,1,mpi_logical,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(is_v_mm,1,mpi_logical,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(is_w_rs_mi,1,mpi_logical,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(is_z_bnd,1,mpi_logical,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(is_ev_qp,1,mpi_logical,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(is_z_qp,1,mpi_logical,0,MPI_COMM_WORLD,
     $    mpierr)
      endif
      if(is_ab_nu_from_tau .eq. .true.) then
        allocate(ab_nu_from_tau(ndim3_tau,0:n_nu_d,2))
      endif
      if(is_tau_from_nu .eq. .true.) then
        allocate(tau_from_nu(0:n_nu,ndim3_tau,2))
      endif
      if(is_ev_bnd .eq. .true.) then
        allocate(ev_bnd(nbasmpw,nbndf,ndim3_k(me3_k+1),nspin_0))
      endif
      if(is_q_qp .eq. .true.) then
        allocate(q_qp(nbndf,nbndf,ndim3_k(me3_k+1),nspin_0))
      endif
      if(is_s_it_r0 .eq. .true.) then
        allocate(s_it_r0(nr_red,n_pbmt,ndim3_kk(me3_k+1)))
      endif
      if(is_v_ii .eq. .true.) then
        allocate(v_ii(nplw_gw,nplw_gw,ndim3_k(me3_k+1)))
      endif
      if(is_v_mm .eq. .true.) then
        allocate(v_mm(nlb2,nlb2,ndim3_k(me3_k+1)))
      endif
      if(is_w_rs_mi .eq. .true.) then
        allocate(w_rs_mi(nr_red,n_pbmt,ndim3_kk(me3_k+1)))
      endif
      if(is_z_bnd .eq. .true.) then
        allocate(z_bnd(nfun,nbndf,ndim3_k(me3_k+1),nspin_0))
      endif
      if(is_ev_qp .eq. .true.) then
        allocate(ev_qp(nbasmpw,nbndf,ndim3_k(me3_k+1),nspin_0))
      endif
      if(is_z_qp .eq. .true.) then
        allocate(z_qp(nfun,nbndf,ndim3_k(me3_k+1),nspin_0))
      endif

      i_len=len_trim(allfile)
c      write(179,*) 'allfile', allfile(1:i_len)      
      allfilet=rstpath(1:len_trim(rstpath))//'/'//
     $  allfile(1:i_len) 
      i_len=len_trim(allfilet)      

      
      if (me3_tau .eq. 0) then
        do ind_k=1,ndim3_k(me3_k+1)
          k0=n3_mpi_k(me3_k+1)+ind_k
          if (k0 .lt. 10) then
            write(k0_ch, '(i1)') k0
          elseif (k0 .lt. 100) then
            write(k0_ch, '(i2)') k0
          elseif (k0 .lt. 1000) then
            write(k0_ch, '(i3)') k0
          elseif (k0 .lt. 10000) then
            write(k0_ch, '(i4)') k0
          endif
          if(is_ev_bnd .eq. .true.) then
            open(3,file=allfilet(1:i_len)//'_ev_bnd_'//trim(k0_ch)
     $        //'_'//trim(ubi)//'.uca',form='unformatted')
            do ispin=1, nspin_0
              read(3)ev_bnd(:,:,ind_k,ispin)
            enddo
            close(3)
          endif
          if(is_q_qp .eq. .true.) then
            open(3,file=allfilet(1:i_len)//'_q_qp_'//trim(k0_ch)
     $        //'_'//trim(ubi)//'.uca',form='unformatted')
            do ispin=1, nspin_0
              read(3)q_qp(:,:,ind_k,ispin)
            enddo
            close(3)
          endif
          if(is_s_it_r0 .eq. .true.) then
            open(3,file=allfilet(1:i_len)//'_s_it_r0_'//trim(k0_ch)
     $        //'_'//trim(ubi)//'.uca',form='unformatted')
            read(3)s_it_r0(:,:,ind_k)
            close(3)
          endif
          if(is_v_ii .eq. .true.) then
            open(3,file=allfilet(1:i_len)//'_v_ii_'//trim(k0_ch)
     $        //'_'//trim(ubi)//'.uca',form='unformatted')
            read(3)v_ii(:,:,ind_k)
            close(3)
          endif
          if(is_v_mm .eq. .true.) then
            open(3,file=allfilet(1:i_len)//'_v_mm_'//trim(k0_ch)
     $        //'_'//trim(ubi)//'.uca',form='unformatted')
            read(3)v_mm(:,:,ind_k)
            close(3)
          endif
          if(is_w_rs_mi .eq. .true.) then
            open(3,file=allfilet(1:i_len)//'_w_rs_mi_'//trim(k0_ch)
     $        //'_'//trim(ubi)//'.uca',form='unformatted')
            read(3)w_rs_mi(:,:,ind_k)
            close(3)
          endif
          if(is_z_bnd .eq. .true.) then
            open(3,file=allfilet(1:i_len)//'_z_bnd_'//trim(k0_ch)
     $        //'_'//trim(ubi)//'.uca',form='unformatted')
            do ispin=1, nspin_0
              read(3)z_bnd(:,:,ind_k,ispin)
            enddo
            close(3)
          endif
          if(is_ev_qp .eq. .true.) then
            open(3,file=allfilet(1:i_len)//'_ev_qp_'//trim(k0_ch)
     $        //'_'//trim(ubi)//'.uca',form='unformatted')
            do ispin=1, nspin_0
              read(3)ev_qp(:,:,ind_k,ispin)
            enddo
            close(3)
          endif
          if(is_z_qp .eq. .true.) then
            open(3,file=allfilet(1:i_len)//'_z_qp_'//trim(k0_ch)
     $        //'_'//trim(ubi)//'.uca',form='unformatted')
            do ispin=1, nspin_0
              read(3)z_qp(:,:,ind_k,ispin)
            enddo
            close(3)
          endif
        enddo
      endif
      i_len=len_trim(allfilet)
      if (me3_k .eq. 0) then
        do ind_tau=1,ndim3_tau
          i_tau=ndim3_tau*me+ind_tau
          if (i_tau .lt. 10) then
            write(tau_ch, '(i1)') i_tau
          elseif (i_tau .lt. 100) then
            write(tau_ch, '(i2)') i_tau
          elseif (i_tau .lt. 1000) then
            write(tau_ch, '(i3)') i_tau
          elseif (i_tau .lt. 10000) then
            write(tau_ch, '(i4)') i_tau
          endif
          if(is_ab_nu_from_tau .eq. .true.) then
            open(3,file=allfilet(1:i_len)//'_ab_nu_from_tau_'
     $        //trim(tau_ch)//'_'//trim(ubi)//'.uca',form='unformatted')
            read(3)ab_nu_from_tau(ind_tau,:,:)
            close(3)
          endif
          if(is_tau_from_nu .eq. .true.) then
            open(3,file=allfilet(1:i_len)//'_tau_from_nu_'//trim(tau_ch)
     $        //'_'//trim(ubi)//'.uca',form='unformatted')
            read(3)tau_from_nu(:,ind_tau,:)
            close(3)
          endif
        enddo
      endif
      i_len=len_trim(allfilet)
      if (me3_k .eq. 0) then
        do ind_nu=1,ndim3_nu
          i_nu=ndim3_nu*me+ind_nu
          if (i_nu .lt. 10) then
            write(nu_ch, '(i1)') i_nu
          elseif (i_nu .lt. 100) then
            write(nu_ch, '(i2)') i_nu
          elseif (i_nu .lt. 1000) then
            write(nu_ch, '(i3)') i_nu
          elseif (i_nu .lt. 10000) then
            write(nu_ch, '(i4)') i_nu
          endif
        enddo
      endif
      i_len=len_trim(allfilet)
      if (me3_tau .eq. 0 .and. me3_k .eq. 0) then
        do ind_k=1,ndim3_k(me3_k+1)
          k0=n3_mpi_k(me3_k+1)+ind_k
          if (k0 .lt. 10) then
            write(k0_ch, '(i1)') k0
          elseif (k0 .lt. 100) then
            write(k0_ch, '(i2)') k0
          elseif (k0 .lt. 1000) then
            write(k0_ch, '(i3)') k0
          elseif (k0 .lt. 10000) then
            write(k0_ch, '(i4)') k0
          endif
          do ind_tau=1,ndim3_tau
            i_tau=ndim3_tau*me+ind_tau
            if (i_tau .lt. 10) then
              write(tau_ch, '(i1)') i_tau
            elseif (i_tau .lt. 100) then
              write(tau_ch, '(i2)') i_tau
            elseif (i_tau .lt. 1000) then
              write(tau_ch, '(i3)') i_tau
            elseif (i_tau .lt. 10000) then
              write(tau_ch, '(i4)') i_tau
            endif
          enddo
        enddo
      endif
      i_len=len_trim(allfilet)
      if (me3_tau .eq. 0 .and. me3_k .eq. 0) then
        do ind_k=1,ndim3_k(me3_k+1)
          k0=n3_mpi_k(me3_k+1)+ind_k
          if (k0 .lt. 10) then
            write(k0_ch, '(i1)') k0
          elseif (k0 .lt. 100) then
            write(k0_ch, '(i2)') k0
          elseif (k0 .lt. 1000) then
            write(k0_ch, '(i3)') k0
          elseif (k0 .lt. 10000) then
            write(k0_ch, '(i4)') k0
          endif
          do ind_nu=1,ndim3_nu
            i_nu=ndim3_tau*me+ind_nu
            if (i_nu .lt. 10) then
              write(nu_ch, '(i1)') i_nu
            elseif (i_nu .lt. 100) then
              write(nu_ch, '(i2)') i_nu
            elseif (i_nu .lt. 1000) then
              write(nu_ch, '(i3)') i_nu
            elseif (i_nu .lt. 10000) then
              write(nu_ch, '(i4)') i_nu
            endif
          enddo
        enddo
      endif

      if(goparr) then
        if(is_ab_nu_from_tau .eq. .true.) then
          call mpi_bcast(ab_nu_from_tau,size(ab_nu_from_tau),
     $      mpi_double_precision,0,comm_tau,mpierr)
        endif
        if(is_tau_from_nu .eq. .true.) then
          call mpi_bcast(tau_from_nu,size(tau_from_nu),
     $      mpi_double_precision,0,comm_tau,mpierr)
        endif
        if(is_ev_bnd .eq. .true.) then
          call mpi_bcast(ev_bnd,size(ev_bnd),mpi_double_complex,0,
     $      comm_k,mpierr)
        endif
        if(is_q_qp .eq. .true.) then
          call mpi_bcast(q_qp,size(q_qp),mpi_double_complex,0,
     $      comm_k,mpierr)
        endif
        if(is_s_it_r0 .eq. .true.) then
          call mpi_bcast(s_it_r0,size(s_it_r0),mpi_double_precision,
     $      0,comm_k,mpierr)
        endif
        if(is_v_ii .eq. .true.) then
          call mpi_bcast(v_ii,size(v_ii),mpi_double_complex,0,
     $      comm_k,mpierr)
        endif
        if(is_v_mm .eq. .true.) then
          call mpi_bcast(v_mm,size(v_mm),mpi_double_complex,0,
     $      comm_k,mpierr)
        endif
        if(is_w_rs_mi .eq. .true.) then
          call mpi_bcast(w_rs_mi,size(w_rs_mi),mpi_double_precision,
     $      0,comm_k,mpierr)
        endif
        if(is_z_bnd .eq. .true.) then
          call mpi_bcast(z_bnd,size(z_bnd),mpi_double_complex,0,
     $      comm_k,mpierr)
        endif
        if(is_ev_qp .eq. .true.) then
          call mpi_bcast(ev_qp,size(ev_qp),mpi_double_complex,0,
     $      comm_k,mpierr)
        endif
        if(is_z_qp .eq. .true.) then
          call mpi_bcast(z_qp,size(z_qp),mpi_double_complex,0,
     $      comm_k,mpierr)
        endif
      endif

      end
