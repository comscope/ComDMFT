      subroutine read_ntau_uca_info
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
      integer :: sz_cg,sz_cgr,sz_indxcg,sz_jcg,sz_jcgr,
     $sz_k_list1,sz_k_list2,sz_ksym_in_npnt
      logical :: is_cg,is_cgr,is_indxcg,is_jcg,is_jcgr,
     $is_k_list,is_ksym_in_npnt
      logical :: is_augm,is_dr,is_gfun,is_gfund,is_h,is_ind_wf,
     $is_indbasa,is_indfun0,is_io_lem,is_iob2,is_iopb,is_lf_isz,
     $is_lfun,is_lfun_pb,is_lfunm,is_li_val,is_lim_pb_mt,is_lm_isz,
     $is_lmb,is_lmpb,is_n_pbmt0,is_nrad,is_ntle,is_ntle_pb,is_r,
     $is_r_atom,is_rl,is_smt,is_txtel,is_bett,is_et_sum_nu,
     $is_ind_nu_int,is_indxcgr,is_tau_mesh,is_w_nu,is_w_nu_d,is_e_bnd,
     $is_e_qp,is_g_sym_0,is_gbs,is_gbs_sym,is_i_kref,is_ifia,
     $is_ig_pair,is_ind_k_star_back,is_ind_r_star,is_ind_r_star_back,
     $is_indgb,is_indplw,is_indpw_gw,is_indpw_r,is_ip,is_ip_k,is_ip_r,
     $is_iplf_gk,is_ir_ref,is_is,is_itshift_r,is_k_group,is_k_star,
     $is_k_sym_0,is_list_sym_k,is_n_bnd,is_nbask,is_nhsort,is_nplw_r,
     $is_nplwgw,is_nrr_red,is_num_r_star,is_num_sym_k,is_pnt,
     $is_r0_pair,is_r_group,is_r_pnt,is_r_pnt_ind,is_shift,is_sovi,
     $is_sovr,is_sym_label,is_tau,is_tshift,is_tshift_r,is_u,is_uj,
     $is_wgt,is_igbs,is_iplf_bk
      logical :: is_ab_nu_from_tau,is_tau_from_nu,is_ev_bnd,
     $is_q_qp,is_s_it_r0,is_v_ii,is_v_mm,is_w_rs_mi,is_z_bnd,is_ev_qp,
     $is_z_qp
      if(maswrk) then
c$$$        write(iun,*) adjustl(trim(rstpath))//'/checkpoint/uca_info.rst'
        open(1,file=adjustl(trim(rstpath))//'/checkpoint/'//
     $    'info.uca',
     $    status='old',form='unformatted')
      read(1)limlb,limlb1,limlpb,max_e_pb,maxb,maxel,maxel_loc,
     $maxlfpb,maxlfun,maxnrad,maxpb,maxpbr,maxntle,maxwf,max_rl,
     $n_pbmt,n_pbmtm,n_pbtot,nfun,nlb2,nrmax,e_c_solid_pw,
     $e_x_solid_pw,allfile,betta_t,chem_pot,chem_pot_qp,complex_ro,
     $eps_pb,interp_nu_d,irel,itest,magn,maxylm,n_omega,n_nu,n_nu_d,
     $n_nu_exa,n_tau,nrel,nrel_core,nspin,nspin_0,nspin_1,ubi,amega,
     $gbas,k_star_max,max_r_star,maxwig,maxwig_j,metal,n_symb,natom,
     $nbasmpw,nbndf,ndiv,ndtau,nkd,nkr,nplw_gw,nplw_rm,npnt,nqdiv,
     $nr_full,nr_red,nrdiv,nrr_max,nrr_max_cor,nsort,numplw,par,r0b0,
     $rb0,rbas,symmorph,v_gg_00,ngroup,maxplw,nplwdif,nbndtm,maxia,
     $maxib,maxic,nplwbas,clight
      endif
      if(goparr) then
      call mpi_bcast(n_tau,1,mpi_integer,0,MPI_COMM_WORLD,
     $mpierr)
      call mpi_bcast(n_nu,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)
      call mpi_bcast(n_nu_asy,1,mpi_integer,0,MPI_COMM_WORLD,
     $mpierr)
      call mpi_bcast(n_nu_exa,1,mpi_integer,0,MPI_COMM_WORLD,
     $mpierr)
      call mpi_bcast(n_omega_exa,1,mpi_integer,0,MPI_COMM_WORLD,
     $mpierr)
      call mpi_bcast(n_omega_asy,1,mpi_integer,0,MPI_COMM_WORLD,
     $mpierr)
      call mpi_bcast(n_omega,1,mpi_integer,0,MPI_COMM_WORLD,
     $mpierr)

      endif
      close(1)
      end
