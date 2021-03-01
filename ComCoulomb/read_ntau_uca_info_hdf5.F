      subroutine read_ntau_uca_info_hdf5
      use manager_mod
      use atom_mod
      use solid_mod
      use parallel_mod
      use units_mod
      use etot_mod
      use wannier90_mod
      use hdf5io_mod
      use hdf5
      implicit none
      include 'mpif.h'
      integer :: i_len,i_lenf,k0,ind_k,ispin
      integer :: i_tau,ind_tau
      character*300 :: tau_ch
      integer :: i_nu,ind_nu
      character*300 :: nu_ch
      character*300 :: k0_ch
      integer(hid_t) :: fid_dat = -1
      integer :: h5err,ierr
      integer :: ng(7)
      integer :: nv(7)
      integer :: itemp
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
c$$$      i_len=len_trim(allfile)
c      if(maswrk) write(iun,*)'ntau before fid_dat',fid_dat
      call finit()
      i_len=len_trim(rstpath)
c$$$      i_lenf=len_trim(allfile)

c$$$      print *, rstpath(1:i_len)//'/info.uca'
      call fopen(rstpath(1:i_len)//'/info.uca',fid_dat)
     
      if(maswrk) write(iun,*) 'ntau after fid_dat',fid_dat
      call fread(MPI_COMM_WORLD,limlb,"limlb",fid_dat)
      call fread(MPI_COMM_WORLD,limlb1,"limlb1",fid_dat)
      call fread(MPI_COMM_WORLD,limlpb,"limlpb",fid_dat)
      call fread(MPI_COMM_WORLD,max_e_pb,"max_e_pb",fid_dat)
      call fread(MPI_COMM_WORLD,maxb,"maxb",fid_dat)
      call fread(MPI_COMM_WORLD,maxel,"maxel",fid_dat)
      call fread(MPI_COMM_WORLD,maxel_loc,"maxel_loc",fid_dat)
      call fread(MPI_COMM_WORLD,maxlfpb,"maxlfpb",fid_dat)
      call fread(MPI_COMM_WORLD,maxlfun,"maxlfun",fid_dat)
      call fread(MPI_COMM_WORLD,maxnrad,"maxnrad",fid_dat)
      call fread(MPI_COMM_WORLD,maxpb,"maxpb",fid_dat)
      call fread(MPI_COMM_WORLD,maxpbr,"maxpbr",fid_dat)
      call fread(MPI_COMM_WORLD,maxntle,"maxntle",fid_dat)
      call fread(MPI_COMM_WORLD,maxwf,"maxwf",fid_dat)
      call fread(MPI_COMM_WORLD,max_rl,"max_rl",fid_dat)
      call fread(MPI_COMM_WORLD,n_pbmt,"n_pbmt",fid_dat)
      call fread(MPI_COMM_WORLD,n_pbmtm,"n_pbmtm",fid_dat)
      call fread(MPI_COMM_WORLD,n_pbtot,"n_pbtot",fid_dat)
      call fread(MPI_COMM_WORLD,nfun,"nfun",fid_dat)
      call fread(MPI_COMM_WORLD,nlb2,"nlb2",fid_dat)
      call fread(MPI_COMM_WORLD,nrmax,"nrmax",fid_dat)
      call fread(MPI_COMM_WORLD,e_c_solid_pw,"e_c_solid_pw",
     $fid_dat)
      call fread(MPI_COMM_WORLD,e_x_solid_pw,"e_x_solid_pw",
     $fid_dat)
      call fread(MPI_COMM_WORLD,allfile,"allfile",fid_dat)
      call fread(MPI_COMM_WORLD,betta_t,"betta_t",fid_dat)
      call fread(MPI_COMM_WORLD,chem_pot,"chem_pot",fid_dat)
      call fread(MPI_COMM_WORLD,chem_pot_qp,"chem_pot_qp",
     $fid_dat)
      call fread(MPI_COMM_WORLD,complex_ro,"complex_ro",fid_dat)
      call fread(MPI_COMM_WORLD,eps_pb,"eps_pb",fid_dat)
      call fread(MPI_COMM_WORLD,interp_nu_d,"interp_nu_d",
     $fid_dat)
      call fread(MPI_COMM_WORLD,irel,"irel",fid_dat)
      call fread(MPI_COMM_WORLD,itest,"itest",fid_dat)
      call fread(MPI_COMM_WORLD,magn,"magn",fid_dat)
      call fread(MPI_COMM_WORLD,maxylm,"maxylm",fid_dat)
      call fread(MPI_COMM_WORLD,n_omega,"n_omega",fid_dat)
      call fread(MPI_COMM_WORLD,n_nu,"n_nu",fid_dat)
      call fread(MPI_COMM_WORLD,n_nu_d,"n_nu_d",fid_dat)
      call fread(MPI_COMM_WORLD,n_nu_exa,"n_nu_exa",fid_dat)
      call fread(MPI_COMM_WORLD,n_tau,"n_tau",fid_dat)
      call fread(MPI_COMM_WORLD,nrel,"nrel",fid_dat)
      call fread(MPI_COMM_WORLD,nrel_core,"nrel_core",fid_dat)
      call fread(MPI_COMM_WORLD,nspin,"nspin",fid_dat)
      call fread(MPI_COMM_WORLD,nspin_0,"nspin_0",fid_dat)
      call fread(MPI_COMM_WORLD,nspin_1,"nspin_1",fid_dat)
      call fread(MPI_COMM_WORLD,ubi,"ubi",fid_dat)
      call fread(MPI_COMM_WORLD,amega,"amega",fid_dat)
      call fread(MPI_COMM_WORLD,k_star_max,"k_star_max",fid_dat)
      call fread(MPI_COMM_WORLD,max_r_star,"max_r_star",fid_dat)
      call fread(MPI_COMM_WORLD,maxwig,"maxwig",fid_dat)
      call fread(MPI_COMM_WORLD,maxwig_j,"maxwig_j",fid_dat)
      call fread(MPI_COMM_WORLD,metal,"metal",fid_dat)
      call fread(MPI_COMM_WORLD,n_symb,"n_symb",fid_dat)
      call fread(MPI_COMM_WORLD,natom,"natom",fid_dat)
      call fread(MPI_COMM_WORLD,nbasmpw,"nbasmpw",fid_dat)
      call fread(MPI_COMM_WORLD,nbndf,"nbndf",fid_dat)
      call fread(MPI_COMM_WORLD,ndtau,"ndtau",fid_dat)
      call fread(MPI_COMM_WORLD,nkd,"nkd",fid_dat)
      call fread(MPI_COMM_WORLD,nkr,"nkr",fid_dat)
      call fread(MPI_COMM_WORLD,nplw_gw,"nplw_gw",fid_dat)
      call fread(MPI_COMM_WORLD,nplw_rm,"nplw_rm",fid_dat)
      call fread(MPI_COMM_WORLD,npnt,"npnt",fid_dat)
      call fread(MPI_COMM_WORLD,nqdiv,"nqdiv",fid_dat)
      call fread(MPI_COMM_WORLD,nr_full,"nr_full",fid_dat)
      call fread(MPI_COMM_WORLD,nr_red,"nr_red",fid_dat)
      call fread(MPI_COMM_WORLD,nrr_max,"nrr_max",fid_dat)
      call fread(MPI_COMM_WORLD,nrr_max_cor,"nrr_max_cor",
     $fid_dat)
      call fread(MPI_COMM_WORLD,nsort,"nsort",fid_dat)
      call fread(MPI_COMM_WORLD,numplw,"numplw",fid_dat)
      call fread(MPI_COMM_WORLD,par,"par",fid_dat)
      call fread(MPI_COMM_WORLD,symmorph,"symmorph",fid_dat)
      call fread(MPI_COMM_WORLD,v_gg_00,"v_gg_00",fid_dat)
      call fread(MPI_COMM_WORLD,ngroup,"ngroup",fid_dat)
      call fread(MPI_COMM_WORLD,maxplw,"maxplw",fid_dat)
      call fread(MPI_COMM_WORLD,nplwdif,"nplwdif",fid_dat)
      call fread(MPI_COMM_WORLD,nbndtm,"nbndtm",fid_dat)
      call fread(MPI_COMM_WORLD,maxia,"maxia",fid_dat)
      call fread(MPI_COMM_WORLD,maxib,"maxib",fid_dat)
      call fread(MPI_COMM_WORLD,maxic,"maxic",fid_dat)
      call fread(MPI_COMM_WORLD,nplwbas,"nplwbas",fid_dat)
      call fread(MPI_COMM_WORLD,clight,"clight",fid_dat)
      if(maswrk) write(iun,*)'read ntau data'
      call fclose(fid_dat)
      end
