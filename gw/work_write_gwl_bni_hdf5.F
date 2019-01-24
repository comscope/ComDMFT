      subroutine write_gwl_bni_hdf5
      use manager_mod
      use atom_mod
      use solid_mod
      use parallel_mod
      use units_mod
      use files_mod
c      use vertex_mod
c      use wannier90_mod
      use etot_mod
c      use models_mod
      use hdf5io_mod
      use hdf5
      implicit none
      include 'mpif.h'
      integer :: i_len,k0,ind_k,ispin
      integer :: i_tau,ind_tau
      character*300 :: tau_ch
      integer :: i_nu,ind_nu
      character*300 :: nu_ch
      character*300 :: k0_ch
c      integer(hid_t) :: fid_dat = -1
      integer :: h5err
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
      i_len=len_trim(allfile)
      write(iun,*)'fid_dat',fid_dat
c       call fopen('gwl.rst',fid_dat)
c       write(iun,*)'fid_dat',fid_dat
      call fwrite(MPI_COMM_WORLD,limlb,"limlb",fid_dat)
      call fwrite(MPI_COMM_WORLD,limlb1,"limlb1",fid_dat)
      call fwrite(MPI_COMM_WORLD,limlpb,"limlpb",fid_dat)
      call fwrite(MPI_COMM_WORLD,max_e_pb,"max_e_pb",fid_dat)
      call fwrite(MPI_COMM_WORLD,maxb,"maxb",fid_dat)
      call fwrite(MPI_COMM_WORLD,maxel,"maxel",fid_dat)
      call fwrite(MPI_COMM_WORLD,maxel_loc,"maxel_loc",fid_dat)
      call fwrite(MPI_COMM_WORLD,maxlfpb,"maxlfpb",fid_dat)
      call fwrite(MPI_COMM_WORLD,maxlfun,"maxlfun",fid_dat)
      call fwrite(MPI_COMM_WORLD,maxnrad,"maxnrad",fid_dat)
      call fwrite(MPI_COMM_WORLD,maxpb,"maxpb",fid_dat)
      call fwrite(MPI_COMM_WORLD,maxpbr,"maxpbr",fid_dat)
      call fwrite(MPI_COMM_WORLD,maxntle,"maxntle",fid_dat)
      call fwrite(MPI_COMM_WORLD,maxwf,"maxwf",fid_dat)
      call fwrite(MPI_COMM_WORLD,max_rl,"max_rl",fid_dat)
      call fwrite(MPI_COMM_WORLD,n_pbmt,"n_pbmt",fid_dat)
      call fwrite(MPI_COMM_WORLD,n_pbmtm,"n_pbmtm",fid_dat)
      call fwrite(MPI_COMM_WORLD,n_pbtot,"n_pbtot",fid_dat)
      call fwrite(MPI_COMM_WORLD,nfun,"nfun",fid_dat)
      call fwrite(MPI_COMM_WORLD,nlb2,"nlb2",fid_dat)
      call fwrite(MPI_COMM_WORLD,nrmax,"nrmax",fid_dat)
      call fwrite(MPI_COMM_WORLD,e_c_solid_pw,"e_c_solid_pw",
     $fid_dat)
      call fwrite(MPI_COMM_WORLD,e_x_solid_pw,"e_x_solid_pw",
     $fid_dat)
      call fwrite(MPI_COMM_WORLD,allfile,"allfile",fid_dat)
      call fwrite(MPI_COMM_WORLD,betta_t,"betta_t",fid_dat)
      call fwrite(MPI_COMM_WORLD,chem_pot,"chem_pot",fid_dat)
      call fwrite(MPI_COMM_WORLD,chem_pot_qp,"chem_pot_qp",
     $fid_dat)
      call fwrite(MPI_COMM_WORLD,complex_ro,"complex_ro",
     $fid_dat)
      call fwrite(MPI_COMM_WORLD,eps_pb,"eps_pb",fid_dat)
      call fwrite(MPI_COMM_WORLD,interp_nu_d,"interp_nu_d",
     $fid_dat)
      call fwrite(MPI_COMM_WORLD,irel,"irel",fid_dat)
      call fwrite(MPI_COMM_WORLD,itest,"itest",fid_dat)
      call fwrite(MPI_COMM_WORLD,magn,"magn",fid_dat)
      call fwrite(MPI_COMM_WORLD,maxylm,"maxylm",fid_dat)
      call fwrite(MPI_COMM_WORLD,n_omega,"n_omega",fid_dat)
      call fwrite(MPI_COMM_WORLD,n_nu,"n_nu",fid_dat)
      call fwrite(MPI_COMM_WORLD,n_nu_d,"n_nu_d",fid_dat)
      call fwrite(MPI_COMM_WORLD,n_nu_exa,"n_nu_exa",fid_dat)
      call fwrite(MPI_COMM_WORLD,n_tau,"n_tau",fid_dat)
      call fwrite(MPI_COMM_WORLD,nrel,"nrel",fid_dat)
      call fwrite(MPI_COMM_WORLD,nrel_core,"nrel_core",fid_dat)
      call fwrite(MPI_COMM_WORLD,nspin,"nspin",fid_dat)
      call fwrite(MPI_COMM_WORLD,nspin_0,"nspin_0",fid_dat)
      call fwrite(MPI_COMM_WORLD,nspin_1,"nspin_1",fid_dat)
      call fwrite(MPI_COMM_WORLD,ubi,"ubi",fid_dat)
      call fwrite(MPI_COMM_WORLD,amega,"amega",fid_dat)
      call fwrite(MPI_COMM_WORLD,k_star_max,"k_star_max",
     $fid_dat)
      call fwrite(MPI_COMM_WORLD,max_r_star,"max_r_star",
     $fid_dat)
      call fwrite(MPI_COMM_WORLD,maxwig,"maxwig",fid_dat)
      call fwrite(MPI_COMM_WORLD,maxwig_j,"maxwig_j",fid_dat)
      call fwrite(MPI_COMM_WORLD,metal,"metal",fid_dat)
      call fwrite(MPI_COMM_WORLD,n_symb,"n_symb",fid_dat)
      call fwrite(MPI_COMM_WORLD,natom,"natom",fid_dat)
      call fwrite(MPI_COMM_WORLD,nbasmpw,"nbasmpw",fid_dat)
      call fwrite(MPI_COMM_WORLD,nbndf,"nbndf",fid_dat)
      call fwrite(MPI_COMM_WORLD,ndtau,"ndtau",fid_dat)
      call fwrite(MPI_COMM_WORLD,nkd,"nkd",fid_dat)
      call fwrite(MPI_COMM_WORLD,nkr,"nkr",fid_dat)
      call fwrite(MPI_COMM_WORLD,nplw_gw,"nplw_gw",fid_dat)
      call fwrite(MPI_COMM_WORLD,nplw_rm,"nplw_rm",fid_dat)
      call fwrite(MPI_COMM_WORLD,npnt,"npnt",fid_dat)
      call fwrite(MPI_COMM_WORLD,nqdiv,"nqdiv",fid_dat)
      call fwrite(MPI_COMM_WORLD,nr_full,"nr_full",fid_dat)
      call fwrite(MPI_COMM_WORLD,nr_red,"nr_red",fid_dat)
      call fwrite(MPI_COMM_WORLD,nrr_max,"nrr_max",fid_dat)
      call fwrite(MPI_COMM_WORLD,nrr_max_cor,"nrr_max_cor",
     $fid_dat)
      call fwrite(MPI_COMM_WORLD,nsort,"nsort",fid_dat)
      call fwrite(MPI_COMM_WORLD,numplw,"numplw",fid_dat)
      call fwrite(MPI_COMM_WORLD,par,"par",fid_dat)
      call fwrite(MPI_COMM_WORLD,symmorph,"symmorph",fid_dat)
      call fwrite(MPI_COMM_WORLD,v_gg_00,"v_gg_00",fid_dat)
      call fwrite(MPI_COMM_WORLD,ngroup,"ngroup",fid_dat)
      call fwrite(MPI_COMM_WORLD,maxplw,"maxplw",fid_dat)
      call fwrite(MPI_COMM_WORLD,nplwdif,"nplwdif",fid_dat)
      call fwrite(MPI_COMM_WORLD,nbndtm,"nbndtm",fid_dat)
      call fwrite(MPI_COMM_WORLD,maxia,"maxia",fid_dat)
      call fwrite(MPI_COMM_WORLD,maxib,"maxib",fid_dat)
      call fwrite(MPI_COMM_WORLD,maxic,"maxic",fid_dat)
      call fwrite(MPI_COMM_WORLD,nplwbas,"nplwbas",fid_dat)
      call fwrite(MPI_COMM_WORLD,clight,"clight",fid_dat)
      nv = 0
      ng(1)=size(gbas,1)
      ng(2)=size(gbas,2)
      call fwrite(MPI_COMM_WORLD,gbas,ng,nv,"gbas",fid_dat)
      nv = 0
      ng(1)=size(ndiv,1)
      call fwrite(MPI_COMM_WORLD,ndiv,ng(1),nv(1),"ndiv",
     $fid_dat)
      nv = 0
      ng(1)=size(nrdiv,1)
      call fwrite(MPI_COMM_WORLD,nrdiv,ng(1),nv(1),"nrdiv",
     $fid_dat)
      nv = 0
      ng(1)=size(r0b0,1)
      ng(2)=size(r0b0,2)
      call fwrite(MPI_COMM_WORLD,r0b0,ng,nv,"r0b0",fid_dat)
      nv = 0
      ng(1)=size(rb0,1)
      ng(2)=size(rb0,2)
      call fwrite(MPI_COMM_WORLD,rb0,ng,nv,"rb0",fid_dat)
      nv = 0
      ng(1)=size(rbas,1)
      ng(2)=size(rbas,2)
      call fwrite(MPI_COMM_WORLD,rbas,ng,nv,"rbas",fid_dat)
      is_augm=allocated(augm)
      is_dr=allocated(dr)
      is_gfun=allocated(gfun)
      is_gfund=allocated(gfund)
      is_h=allocated(h)
      is_ind_wf=allocated(ind_wf)
      is_indbasa=allocated(indbasa)
      is_indfun0=allocated(indfun0)
      is_io_lem=allocated(io_lem)
      is_iob2=allocated(iob2)
      is_iopb=allocated(iopb)
      is_lf_isz=allocated(lf_isz)
      is_lfun=allocated(lfun)
      is_lfun_pb=allocated(lfun_pb)
      is_lfunm=allocated(lfunm)
      is_li_val=allocated(li_val)
      is_lim_pb_mt=allocated(lim_pb_mt)
      is_lm_isz=allocated(lm_isz)
      is_lmb=allocated(lmb)
      is_lmpb=allocated(lmpb)
      is_n_pbmt0=allocated(n_pbmt0)
      is_nrad=allocated(nrad)
      is_ntle=allocated(ntle)
      is_ntle_pb=allocated(ntle_pb)
      is_r=allocated(r)
      is_r_atom=allocated(r_atom)
      is_rl=allocated(rl)
      is_smt=allocated(smt)
      is_txtel=allocated(txtel)
      is_bett=allocated(bett)
      is_et_sum_nu=allocated(et_sum_nu)
      is_ind_nu_int=allocated(ind_nu_int)
      is_indxcgr=allocated(indxcgr)
      is_tau_mesh=allocated(tau_mesh)
      is_w_nu=allocated(w_nu)
      is_w_nu_d=allocated(w_nu_d)
      is_e_bnd=allocated(e_bnd)
      is_e_qp=allocated(e_qp)
      is_g_sym_0=allocated(g_sym_0)
      is_gbs=allocated(gbs)
      is_gbs_sym=allocated(gbs_sym)
      is_i_kref=allocated(i_kref)
      is_ifia=allocated(ifia)
      is_ig_pair=allocated(ig_pair)
      is_ind_k_star_back=allocated(ind_k_star_back)
      is_ind_r_star=allocated(ind_r_star)
      is_ind_r_star_back=allocated(ind_r_star_back)
      is_indgb=allocated(indgb)
      is_indplw=allocated(indplw)
      is_indpw_gw=allocated(indpw_gw)
      is_indpw_r=allocated(indpw_r)
      is_ip=allocated(ip)
      is_ip_k=allocated(ip_k)
      is_ip_r=allocated(ip_r)
      is_iplf_gk=allocated(iplf_gk)
      is_ir_ref=allocated(ir_ref)
      is_is=allocated(is)
      is_itshift_r=allocated(itshift_r)
      is_k_group=allocated(k_group)
      is_k_star=allocated(k_star)
      is_k_sym_0=allocated(k_sym_0)
      is_list_sym_k=allocated(list_sym_k)
      is_n_bnd=allocated(n_bnd)
      is_nbask=allocated(nbask)
      is_nhsort=allocated(nhsort)
      is_nplw_r=allocated(nplw_r)
      is_nplwgw=allocated(nplwgw)
      is_nrr_red=allocated(nrr_red)
      is_num_r_star=allocated(num_r_star)
      is_num_sym_k=allocated(num_sym_k)
      is_pnt=allocated(pnt)
      is_r0_pair=allocated(r0_pair)
      is_r_group=allocated(r_group)
      is_r_pnt=allocated(r_pnt)
      is_r_pnt_ind=allocated(r_pnt_ind)
      is_shift=allocated(shift)
      is_sovi=allocated(sovi)
      is_sovr=allocated(sovr)
      is_sym_label=allocated(sym_label)
      is_tau=allocated(tau)
      is_tshift=allocated(tshift)
      is_tshift_r=allocated(tshift_r)
      is_u=allocated(u)
      is_uj=allocated(uj)
      is_wgt=allocated(wgt)
      is_igbs=allocated(igbs)
      is_iplf_bk=allocated(iplf_bk)
      call fwrite(MPI_COMM_WORLD,is_augm,"is_augm",fid_dat)
      call fwrite(MPI_COMM_WORLD,is_dr,"is_dr",fid_dat)
      call fwrite(MPI_COMM_WORLD,is_gfun,"is_gfun",fid_dat)
      call fwrite(MPI_COMM_WORLD,is_gfund,"is_gfund",fid_dat)
      call fwrite(MPI_COMM_WORLD,is_h,"is_h",fid_dat)
      call fwrite(MPI_COMM_WORLD,is_ind_wf,"is_ind_wf",fid_dat)
      call fwrite(MPI_COMM_WORLD,is_indbasa,"is_indbasa",
     $fid_dat)
      call fwrite(MPI_COMM_WORLD,is_indfun0,"is_indfun0",
     $fid_dat)
      call fwrite(MPI_COMM_WORLD,is_io_lem,"is_io_lem",fid_dat)
      call fwrite(MPI_COMM_WORLD,is_iob2,"is_iob2",fid_dat)
      call fwrite(MPI_COMM_WORLD,is_iopb,"is_iopb",fid_dat)
      call fwrite(MPI_COMM_WORLD,is_lf_isz,"is_lf_isz",fid_dat)
      call fwrite(MPI_COMM_WORLD,is_lfun,"is_lfun",fid_dat)
      call fwrite(MPI_COMM_WORLD,is_lfun_pb,"is_lfun_pb",
     $fid_dat)
      call fwrite(MPI_COMM_WORLD,is_lfunm,"is_lfunm",fid_dat)
      call fwrite(MPI_COMM_WORLD,is_li_val,"is_li_val",fid_dat)
      call fwrite(MPI_COMM_WORLD,is_lim_pb_mt,"is_lim_pb_mt",
     $fid_dat)
      call fwrite(MPI_COMM_WORLD,is_lm_isz,"is_lm_isz",fid_dat)
      call fwrite(MPI_COMM_WORLD,is_lmb,"is_lmb",fid_dat)
      call fwrite(MPI_COMM_WORLD,is_lmpb,"is_lmpb",fid_dat)
      call fwrite(MPI_COMM_WORLD,is_n_pbmt0,"is_n_pbmt0",
     $fid_dat)
      call fwrite(MPI_COMM_WORLD,is_nrad,"is_nrad",fid_dat)
      call fwrite(MPI_COMM_WORLD,is_ntle,"is_ntle",fid_dat)
      call fwrite(MPI_COMM_WORLD,is_ntle_pb,"is_ntle_pb",
     $fid_dat)
      call fwrite(MPI_COMM_WORLD,is_r,"is_r",fid_dat)
      call fwrite(MPI_COMM_WORLD,is_r_atom,"is_r_atom",fid_dat)
      call fwrite(MPI_COMM_WORLD,is_rl,"is_rl",fid_dat)
      call fwrite(MPI_COMM_WORLD,is_smt,"is_smt",fid_dat)
      call fwrite(MPI_COMM_WORLD,is_txtel,"is_txtel",fid_dat)
      call fwrite(MPI_COMM_WORLD,is_bett,"is_bett",fid_dat)
      call fwrite(MPI_COMM_WORLD,is_et_sum_nu,"is_et_sum_nu",
     $fid_dat)
      call fwrite(MPI_COMM_WORLD,is_ind_nu_int,"is_ind_nu_int",
     $fid_dat)
      call fwrite(MPI_COMM_WORLD,is_indxcgr,"is_indxcgr",
     $fid_dat)
      call fwrite(MPI_COMM_WORLD,is_tau_mesh,"is_tau_mesh",
     $fid_dat)
      call fwrite(MPI_COMM_WORLD,is_w_nu,"is_w_nu",fid_dat)
      call fwrite(MPI_COMM_WORLD,is_w_nu_d,"is_w_nu_d",fid_dat)
      call fwrite(MPI_COMM_WORLD,is_e_bnd,"is_e_bnd",fid_dat)
      call fwrite(MPI_COMM_WORLD,is_e_qp,"is_e_qp",fid_dat)
      call fwrite(MPI_COMM_WORLD,is_g_sym_0,"is_g_sym_0",
     $fid_dat)
      call fwrite(MPI_COMM_WORLD,is_gbs,"is_gbs",fid_dat)
      call fwrite(MPI_COMM_WORLD,is_gbs_sym,"is_gbs_sym",
     $fid_dat)
      call fwrite(MPI_COMM_WORLD,is_i_kref,"is_i_kref",fid_dat)
      call fwrite(MPI_COMM_WORLD,is_ifia,"is_ifia",fid_dat)
      call fwrite(MPI_COMM_WORLD,is_ig_pair,"is_ig_pair",
     $fid_dat)
      call fwrite(MPI_COMM_WORLD,is_ind_k_star_back,
     $"is_ind_k_star_back",fid_dat)
      call fwrite(MPI_COMM_WORLD,is_ind_r_star,"is_ind_r_star",
     $fid_dat)
      call fwrite(MPI_COMM_WORLD,is_ind_r_star_back,
     $"is_ind_r_star_back",fid_dat)
      call fwrite(MPI_COMM_WORLD,is_indgb,"is_indgb",fid_dat)
      call fwrite(MPI_COMM_WORLD,is_indplw,"is_indplw",fid_dat)
      call fwrite(MPI_COMM_WORLD,is_indpw_gw,"is_indpw_gw",
     $fid_dat)
      call fwrite(MPI_COMM_WORLD,is_indpw_r,"is_indpw_r",
     $fid_dat)
      call fwrite(MPI_COMM_WORLD,is_ip,"is_ip",fid_dat)
      call fwrite(MPI_COMM_WORLD,is_ip_k,"is_ip_k",fid_dat)
      call fwrite(MPI_COMM_WORLD,is_ip_r,"is_ip_r",fid_dat)
      call fwrite(MPI_COMM_WORLD,is_iplf_gk,"is_iplf_gk",
     $fid_dat)
      call fwrite(MPI_COMM_WORLD,is_ir_ref,"is_ir_ref",fid_dat)
      call fwrite(MPI_COMM_WORLD,is_is,"is_is",fid_dat)
      call fwrite(MPI_COMM_WORLD,is_itshift_r,"is_itshift_r",
     $fid_dat)
      call fwrite(MPI_COMM_WORLD,is_k_group,"is_k_group",
     $fid_dat)
      call fwrite(MPI_COMM_WORLD,is_k_star,"is_k_star",fid_dat)
      call fwrite(MPI_COMM_WORLD,is_k_sym_0,"is_k_sym_0",
     $fid_dat)
      call fwrite(MPI_COMM_WORLD,is_list_sym_k,"is_list_sym_k",
     $fid_dat)
      call fwrite(MPI_COMM_WORLD,is_n_bnd,"is_n_bnd",fid_dat)
      call fwrite(MPI_COMM_WORLD,is_nbask,"is_nbask",fid_dat)
      call fwrite(MPI_COMM_WORLD,is_nhsort,"is_nhsort",fid_dat)
      call fwrite(MPI_COMM_WORLD,is_nplw_r,"is_nplw_r",fid_dat)
      call fwrite(MPI_COMM_WORLD,is_nplwgw,"is_nplwgw",fid_dat)
      call fwrite(MPI_COMM_WORLD,is_nrr_red,"is_nrr_red",
     $fid_dat)
      call fwrite(MPI_COMM_WORLD,is_num_r_star,"is_num_r_star",
     $fid_dat)
      call fwrite(MPI_COMM_WORLD,is_num_sym_k,"is_num_sym_k",
     $fid_dat)
      call fwrite(MPI_COMM_WORLD,is_pnt,"is_pnt",fid_dat)
      call fwrite(MPI_COMM_WORLD,is_r0_pair,"is_r0_pair",
     $fid_dat)
      call fwrite(MPI_COMM_WORLD,is_r_group,"is_r_group",
     $fid_dat)
      call fwrite(MPI_COMM_WORLD,is_r_pnt,"is_r_pnt",fid_dat)
      call fwrite(MPI_COMM_WORLD,is_r_pnt_ind,"is_r_pnt_ind",
     $fid_dat)
      call fwrite(MPI_COMM_WORLD,is_shift,"is_shift",fid_dat)
      call fwrite(MPI_COMM_WORLD,is_sovi,"is_sovi",fid_dat)
      call fwrite(MPI_COMM_WORLD,is_sovr,"is_sovr",fid_dat)
      call fwrite(MPI_COMM_WORLD,is_sym_label,"is_sym_label",
     $fid_dat)
      call fwrite(MPI_COMM_WORLD,is_tau,"is_tau",fid_dat)
      call fwrite(MPI_COMM_WORLD,is_tshift,"is_tshift",fid_dat)
      call fwrite(MPI_COMM_WORLD,is_tshift_r,"is_tshift_r",
     $fid_dat)
      call fwrite(MPI_COMM_WORLD,is_u,"is_u",fid_dat)
      call fwrite(MPI_COMM_WORLD,is_uj,"is_uj",fid_dat)
      call fwrite(MPI_COMM_WORLD,is_wgt,"is_wgt",fid_dat)
      call fwrite(MPI_COMM_WORLD,is_igbs,"is_igbs",fid_dat)
      call fwrite(MPI_COMM_WORLD,is_iplf_bk,"is_iplf_bk",
     $fid_dat)
      if(is_augm .eq. .true.) then
      nv = 0
      ng(1)=size(augm,1)
      ng(2)=size(augm,2)
      ng(3)=size(augm,3)
      call fwrite(MPI_COMM_WORLD,augm,ng,nv,"augm",fid_dat)
      endif
      if(is_dr .eq. .true.) then
      nv = 0
      ng(1)=size(dr,1)
      ng(2)=size(dr,2)
      call fwrite(MPI_COMM_WORLD,dr,ng,nv,"dr",fid_dat)
      endif
      if(is_gfun .eq. .true.) then
      nv = 0
      ng(1)=size(gfun,1)
      ng(2)=size(gfun,2)
      call fwrite(MPI_COMM_WORLD,gfun,ng,nv,"gfun",fid_dat)
      endif
      if(is_gfund .eq. .true.) then
      nv = 0
      ng(1)=size(gfund,1)
      ng(2)=size(gfund,2)
      call fwrite(MPI_COMM_WORLD,gfund,ng,nv,"gfund",fid_dat)
      endif
      if(is_h .eq. .true.) then
      nv = 0
      ng(1)=size(h,1)
      call fwrite(MPI_COMM_WORLD,h,ng(1),nv(1),"h",fid_dat)
      endif
      if(is_ind_wf .eq. .true.) then
      nv = 0
      ng(1)=size(ind_wf,1)
      ng(2)=size(ind_wf,2)
      call fwrite(MPI_COMM_WORLD,ind_wf,ng,nv,"ind_wf",fid_dat)
      endif
      if(is_indbasa .eq. .true.) then
      nv = 0
      ng(1)=size(indbasa,1)
      ng(2)=size(indbasa,2)
      ng(3)=size(indbasa,3)
      ng(4)=size(indbasa,4)
      call fwrite(MPI_COMM_WORLD,indbasa,ng,nv,"indbasa",
     $fid_dat)
      endif
      if(is_indfun0 .eq. .true.) then
      nv = 0
      ng(1)=size(indfun0,1)
      ng(2)=size(indfun0,2)
      ng(3)=size(indfun0,3)
      ng(4)=size(indfun0,4)
      call fwrite(MPI_COMM_WORLD,indfun0,ng,nv,"indfun0",
     $fid_dat)
      endif
      if(is_io_lem .eq. .true.) then
      nv = 0
      ng(1)=size(io_lem,1)
      call fwrite(MPI_COMM_WORLD,io_lem,ng(1),nv(1),"io_lem",
     $fid_dat)
      endif
      if(is_iob2 .eq. .true.) then
      nv = 0
      ng(1)=size(iob2,1)
      call fwrite(MPI_COMM_WORLD,iob2,ng(1),nv(1),"iob2",
     $fid_dat)
      endif
      if(is_iopb .eq. .true.) then
      nv = 0
      ng(1)=size(iopb,1)
      call fwrite(MPI_COMM_WORLD,iopb,ng(1),nv(1),"iopb",
     $fid_dat)
      endif
      if(is_lf_isz .eq. .true.) then
      nv = 0
      ng(1)=size(lf_isz,1)
      ng(2)=size(lf_isz,2)
      call fwrite(MPI_COMM_WORLD,lf_isz,ng,nv,"lf_isz",fid_dat)
      endif
      if(is_lfun .eq. .true.) then
      nv = 0
      ng(1)=size(lfun,1)
      call fwrite(MPI_COMM_WORLD,lfun,ng(1),nv(1),"lfun",
     $fid_dat)
      endif
      if(is_lfun_pb .eq. .true.) then
      nv = 0
      ng(1)=size(lfun_pb,1)
      call fwrite(MPI_COMM_WORLD,lfun_pb,ng(1),nv(1),"lfun_pb",
     $fid_dat)
      endif
      if(is_lfunm .eq. .true.) then
      nv = 0
      ng(1)=size(lfunm,1)
      call fwrite(MPI_COMM_WORLD,lfunm,ng(1),nv(1),"lfunm",
     $fid_dat)
      endif
      if(is_li_val .eq. .true.) then
      nv = 0
      ng(1)=size(li_val,1)
      ng(2)=size(li_val,2)
      call fwrite(MPI_COMM_WORLD,li_val,ng,nv,"li_val",fid_dat)
      endif
      if(is_lim_pb_mt .eq. .true.) then
      nv = 0
      ng(1)=size(lim_pb_mt,1)
      ng(2)=size(lim_pb_mt,2)
      call fwrite(MPI_COMM_WORLD,lim_pb_mt,ng,nv,"lim_pb_mt",
     $fid_dat)
      endif
      if(is_lm_isz .eq. .true.) then
      nv = 0
      ng(1)=size(lm_isz,1)
      ng(2)=size(lm_isz,2)
      call fwrite(MPI_COMM_WORLD,lm_isz,ng,nv,"lm_isz",fid_dat)
      endif
      if(is_lmb .eq. .true.) then
      nv = 0
      ng(1)=size(lmb,1)
      call fwrite(MPI_COMM_WORLD,lmb,ng(1),nv(1),"lmb",fid_dat)
      endif
      if(is_lmpb .eq. .true.) then
      nv = 0
      ng(1)=size(lmpb,1)
      call fwrite(MPI_COMM_WORLD,lmpb,ng(1),nv(1),"lmpb",
     $fid_dat)
      endif
      if(is_n_pbmt0 .eq. .true.) then
      nv = 0
      ng(1)=size(n_pbmt0,1)
      call fwrite(MPI_COMM_WORLD,n_pbmt0,ng(1),nv(1),"n_pbmt0",
     $fid_dat)
      endif
      if(is_nrad .eq. .true.) then
      nv = 0
      ng(1)=size(nrad,1)
      call fwrite(MPI_COMM_WORLD,nrad,ng(1),nv(1),"nrad",
     $fid_dat)
      endif
      if(is_ntle .eq. .true.) then
      nv = 0
      ng(1)=size(ntle,1)
      ng(2)=size(ntle,2)
      call fwrite(MPI_COMM_WORLD,ntle,ng,nv,"ntle",fid_dat)
      endif
      if(is_ntle_pb .eq. .true.) then
      nv = 0
      ng(1)=size(ntle_pb,1)
      ng(2)=size(ntle_pb,2)
      call fwrite(MPI_COMM_WORLD,ntle_pb,ng,nv,"ntle_pb",
     $fid_dat)
      endif
      if(is_r .eq. .true.) then
      nv = 0
      ng(1)=size(r,1)
      ng(2)=size(r,2)
      call fwrite(MPI_COMM_WORLD,r,ng,nv,"r",fid_dat)
      endif
      if(is_r_atom .eq. .true.) then
      nv = 0
      ng(1)=size(r_atom,1)
      ng(2)=size(r_atom,2)
      call fwrite(MPI_COMM_WORLD,r_atom,ng,nv,"r_atom",fid_dat)
      endif
      if(is_rl .eq. .true.) then
      nv = 0
      ng(1)=size(rl,1)
      ng(2)=size(rl,2)
      ng(3)=size(rl,3)
      call fwrite(MPI_COMM_WORLD,rl,ng,nv,"rl",fid_dat)
      endif
      if(is_smt .eq. .true.) then
      nv = 0
      ng(1)=size(smt,1)
      call fwrite(MPI_COMM_WORLD,smt,ng(1),nv(1),"smt",fid_dat)
      endif
      if(is_txtel .eq. .true.) then
      nv = 0
      ng(1)=size(txtel,1)
      call fwrite(MPI_COMM_WORLD,txtel,ng(1),nv(1),"txtel",
     $fid_dat)
      endif
      if(is_bett .eq. .true.) then
      nv = 0
      ng(1)=size(bett,1)
      call fwrite(MPI_COMM_WORLD,bett,ng(1),nv(1),"bett",
     $fid_dat)
      endif
      if(is_et_sum_nu .eq. .true.) then
      nv = 0
      ng(1)=size(et_sum_nu,1)
      call fwrite(MPI_COMM_WORLD,et_sum_nu,ng(1),nv(1),
     $"et_sum_nu",fid_dat)
      endif
      if(is_ind_nu_int .eq. .true.) then
      nv = 0
      ng(1)=size(ind_nu_int,1)
      call fwrite(MPI_COMM_WORLD,ind_nu_int,ng(1),nv(1),
     $"ind_nu_int",fid_dat)
      endif
      if(is_indxcgr .eq. .true.) then
      nv = 0
      ng(1)=size(indxcgr,1)
      ng(2)=size(indxcgr,2)
      ng(3)=size(indxcgr,3)
      ng(4)=size(indxcgr,4)
      call fwrite(MPI_COMM_WORLD,indxcgr,ng,nv,"indxcgr",
     $fid_dat)
      endif
      if(is_tau_mesh .eq. .true.) then
      nv = 0
      ng(1)=size(tau_mesh,1)
      call fwrite(MPI_COMM_WORLD,tau_mesh,ng(1),nv(1),
     $"tau_mesh",fid_dat)
      endif
      if(is_w_nu .eq. .true.) then
      nv = 0
      ng(1)=size(w_nu,1)
      call fwrite(MPI_COMM_WORLD,w_nu,ng(1),nv(1),"w_nu",
     $fid_dat)
      endif
      if(is_w_nu_d .eq. .true.) then
      nv = 0
      ng(1)=size(w_nu_d,1)
      call fwrite(MPI_COMM_WORLD,w_nu_d,ng(1),nv(1),"w_nu_d",
     $fid_dat)
      endif
      if(is_e_bnd .eq. .true.) then
      nv = 0
      ng(1)=size(e_bnd,1)
      ng(2)=size(e_bnd,2)
      ng(3)=size(e_bnd,3)
      call fwrite(MPI_COMM_WORLD,e_bnd,ng,nv,"e_bnd",fid_dat)
      endif
      if(is_e_qp .eq. .true.) then
      nv = 0
      ng(1)=size(e_qp,1)
      ng(2)=size(e_qp,2)
      ng(3)=size(e_qp,3)
      call fwrite(MPI_COMM_WORLD,e_qp,ng,nv,"e_qp",fid_dat)
      endif
      if(is_g_sym_0 .eq. .true.) then
      nv = 0
      ng(1)=size(g_sym_0,1)
      ng(2)=size(g_sym_0,2)
      call fwrite(MPI_COMM_WORLD,g_sym_0,ng,nv,"g_sym_0",
     $fid_dat)
      endif
      if(is_gbs .eq. .true.) then
      nv = 0
      ng(1)=size(gbs,1)
      ng(2)=size(gbs,2)
      call fwrite(MPI_COMM_WORLD,gbs,ng,nv,"gbs",fid_dat)
      endif
      if(is_gbs_sym .eq. .true.) then
      nv = 0
      ng(1)=size(gbs_sym,1)
      ng(2)=size(gbs_sym,2)
      call fwrite(MPI_COMM_WORLD,gbs_sym,ng,nv,"gbs_sym",
     $fid_dat)
      endif
      if(is_i_kref .eq. .true.) then
      nv = 0
      ng(1)=size(i_kref,1)
      call fwrite(MPI_COMM_WORLD,i_kref,ng(1),nv(1),"i_kref",
     $fid_dat)
      endif
      if(is_ifia .eq. .true.) then
      nv = 0
      ng(1)=size(ifia,1)
      call fwrite(MPI_COMM_WORLD,ifia,ng(1),nv(1),"ifia",
     $fid_dat)
      endif
      if(is_ig_pair .eq. .true.) then
      nv = 0
      ng(1)=size(ig_pair,1)
      ng(2)=size(ig_pair,2)
      ng(3)=size(ig_pair,3)
      call fwrite(MPI_COMM_WORLD,ig_pair,ng,nv,"ig_pair",
     $fid_dat)
      endif
      if(is_ind_k_star_back .eq. .true.) then
      nv = 0
      ng(1)=size(ind_k_star_back,1)
      call fwrite(MPI_COMM_WORLD,ind_k_star_back,ng(1),nv(1),
     $"ind_k_star_back",fid_dat)
      endif
      if(is_ind_r_star .eq. .true.) then
      nv = 0
      ng(1)=size(ind_r_star,1)
      ng(2)=size(ind_r_star,2)
      call fwrite(MPI_COMM_WORLD,ind_r_star,ng,nv,"ind_r_star",
     $fid_dat)
      endif
      if(is_ind_r_star_back .eq. .true.) then
      nv = 0
      ng(1)=size(ind_r_star_back,1)
      call fwrite(MPI_COMM_WORLD,ind_r_star_back,ng(1),nv(1),
     $"ind_r_star_back",fid_dat)
      endif
      if(is_indgb .eq. .true.) then
      nv = 0
      ng(1)=size(indgb,1)
      ng(2)=size(indgb,2)
      call fwrite(MPI_COMM_WORLD,indgb,ng,nv,"indgb",fid_dat)
      endif
      if(is_indplw .eq. .true.) then
      nv = 0
      ng(1)=size(indplw,1)
      ng(2)=size(indplw,2)
      ng(3)=size(indplw,3)
      call fwrite(MPI_COMM_WORLD,indplw,ng,nv,"indplw",fid_dat)
      endif
      if(is_indpw_gw .eq. .true.) then
      nv = 0
      ng(1)=size(indpw_gw,1)
      ng(2)=size(indpw_gw,2)
      call fwrite(MPI_COMM_WORLD,indpw_gw,ng,nv,"indpw_gw",
     $fid_dat)
      endif
      if(is_indpw_r .eq. .true.) then
      nv = 0
      ng(1)=size(indpw_r,1)
      ng(2)=size(indpw_r,2)
      call fwrite(MPI_COMM_WORLD,indpw_r,ng,nv,"indpw_r",
     $fid_dat)
      endif
      if(is_ip .eq. .true.) then
      nv = 0
      ng(1)=size(ip,1)
      ng(2)=size(ip,2)
      call fwrite(MPI_COMM_WORLD,ip,ng,nv,"ip",fid_dat)
      endif
      if(is_ip_k .eq. .true.) then
      nv = 0
      ng(1)=size(ip_k,1)
      ng(2)=size(ip_k,2)
      call fwrite(MPI_COMM_WORLD,ip_k,ng,nv,"ip_k",fid_dat)
      endif
      if(is_ip_r .eq. .true.) then
      nv = 0
      ng(1)=size(ip_r,1)
      ng(2)=size(ip_r,2)
      call fwrite(MPI_COMM_WORLD,ip_r,ng,nv,"ip_r",fid_dat)
      endif
      if(is_iplf_gk .eq. .true.) then
      nv = 0
      ng(1)=size(iplf_gk,1)
      ng(2)=size(iplf_gk,2)
      call fwrite(MPI_COMM_WORLD,iplf_gk,ng,nv,"iplf_gk",
     $fid_dat)
      endif
      if(is_ir_ref .eq. .true.) then
      nv = 0
      ng(1)=size(ir_ref,1)
      call fwrite(MPI_COMM_WORLD,ir_ref,ng(1),nv(1),"ir_ref",
     $fid_dat)
      endif
      if(is_is .eq. .true.) then
      nv = 0
      ng(1)=size(is,1)
      call fwrite(MPI_COMM_WORLD,is,ng(1),nv(1),"is",fid_dat)
      endif
      if(is_itshift_r .eq. .true.) then
      nv = 0
      ng(1)=size(itshift_r,1)
      ng(2)=size(itshift_r,2)
      ng(3)=size(itshift_r,3)
      call fwrite(MPI_COMM_WORLD,itshift_r,ng,nv,"itshift_r",
     $fid_dat)
      endif
      if(is_k_group .eq. .true.) then
      nv = 0
      ng(1)=size(k_group,1)
      call fwrite(MPI_COMM_WORLD,k_group,ng(1),nv(1),"k_group",
     $fid_dat)
      endif
      if(is_k_star .eq. .true.) then
      nv = 0
      ng(1)=size(k_star,1)
      call fwrite(MPI_COMM_WORLD,k_star,ng(1),nv(1),"k_star",
     $fid_dat)
      endif
      if(is_k_sym_0 .eq. .true.) then
      nv = 0
      ng(1)=size(k_sym_0,1)
      ng(2)=size(k_sym_0,2)
      call fwrite(MPI_COMM_WORLD,k_sym_0,ng,nv,"k_sym_0",
     $fid_dat)
      endif
      if(is_list_sym_k .eq. .true.) then
      nv = 0
      ng(1)=size(list_sym_k,1)
      ng(2)=size(list_sym_k,2)
      call fwrite(MPI_COMM_WORLD,list_sym_k,ng,nv,"list_sym_k",
     $fid_dat)
      endif
      if(is_n_bnd .eq. .true.) then
      nv = 0
      ng(1)=size(n_bnd,1)
      ng(2)=size(n_bnd,2)
      call fwrite(MPI_COMM_WORLD,n_bnd,ng,nv,"n_bnd",fid_dat)
      endif
      if(is_nbask .eq. .true.) then
      nv = 0
      ng(1)=size(nbask,1)
      call fwrite(MPI_COMM_WORLD,nbask,ng(1),nv(1),"nbask",
     $fid_dat)
      endif
      if(is_nhsort .eq. .true.) then
      nv = 0
      ng(1)=size(nhsort,1)
      call fwrite(MPI_COMM_WORLD,nhsort,ng(1),nv(1),"nhsort",
     $fid_dat)
      endif
      if(is_nplw_r .eq. .true.) then
      nv = 0
      ng(1)=size(nplw_r,1)
      call fwrite(MPI_COMM_WORLD,nplw_r,ng(1),nv(1),"nplw_r",
     $fid_dat)
      endif
      if(is_nplwgw .eq. .true.) then
      nv = 0
      ng(1)=size(nplwgw,1)
      call fwrite(MPI_COMM_WORLD,nplwgw,ng(1),nv(1),"nplwgw",
     $fid_dat)
      endif
      if(is_nrr_red .eq. .true.) then
      nv = 0
      ng(1)=size(nrr_red,1)
      ng(2)=size(nrr_red,2)
      call fwrite(MPI_COMM_WORLD,nrr_red,ng,nv,"nrr_red",
     $fid_dat)
      endif
      if(is_num_r_star .eq. .true.) then
      nv = 0
      ng(1)=size(num_r_star,1)
      call fwrite(MPI_COMM_WORLD,num_r_star,ng(1),nv(1),
     $"num_r_star",fid_dat)
      endif
      if(is_num_sym_k .eq. .true.) then
      nv = 0
      ng(1)=size(num_sym_k,1)
      call fwrite(MPI_COMM_WORLD,num_sym_k,ng(1),nv(1),
     $"num_sym_k",fid_dat)
      endif
      if(is_pnt .eq. .true.) then
      nv = 0
      ng(1)=size(pnt,1)
      ng(2)=size(pnt,2)
      call fwrite(MPI_COMM_WORLD,pnt,ng,nv,"pnt",fid_dat)
      endif
      if(is_r0_pair .eq. .true.) then
      nv = 0
      ng(1)=size(r0_pair,1)
      ng(2)=size(r0_pair,2)
      ng(3)=size(r0_pair,3)
      call fwrite(MPI_COMM_WORLD,r0_pair,ng,nv,"r0_pair",
     $fid_dat)
      endif
      if(is_r_group .eq. .true.) then
      nv = 0
      ng(1)=size(r_group,1)
      call fwrite(MPI_COMM_WORLD,r_group,ng(1),nv(1),"r_group",
     $fid_dat)
      endif
      if(is_r_pnt .eq. .true.) then
      nv = 0
      ng(1)=size(r_pnt,1)
      ng(2)=size(r_pnt,2)
      call fwrite(MPI_COMM_WORLD,r_pnt,ng,nv,"r_pnt",fid_dat)
      endif
      if(is_r_pnt_ind .eq. .true.) then
      nv = 0
      ng(1)=size(r_pnt_ind,1)
      call fwrite(MPI_COMM_WORLD,r_pnt_ind,ng(1),nv(1),
     $"r_pnt_ind",fid_dat)
      endif
      if(is_shift .eq. .true.) then
      nv = 0
      ng(1)=size(shift,1)
      ng(2)=size(shift,2)
      call fwrite(MPI_COMM_WORLD,shift,ng,nv,"shift",fid_dat)
      endif
      if(is_sovi .eq. .true.) then
      nv = 0
      ng(1)=size(sovi,1)
      call fwrite(MPI_COMM_WORLD,sovi,ng(1),nv(1),"sovi",
     $fid_dat)
      endif
      if(is_sovr .eq. .true.) then
      nv = 0
      ng(1)=size(sovr,1)
      call fwrite(MPI_COMM_WORLD,sovr,ng(1),nv(1),"sovr",
     $fid_dat)
      endif
      if(is_sym_label .eq. .true.) then
      nv = 0
      ng(1)=size(sym_label,1)
      call fwrite(MPI_COMM_WORLD,sym_label,ng(1),nv(1),
     $"sym_label",fid_dat)
      endif
      if(is_tau .eq. .true.) then
      nv = 0
      ng(1)=size(tau,1)
      ng(2)=size(tau,2)
      call fwrite(MPI_COMM_WORLD,tau,ng,nv,"tau",fid_dat)
      endif
      if(is_tshift .eq. .true.) then
      nv = 0
      ng(1)=size(tshift,1)
      ng(2)=size(tshift,2)
      ng(3)=size(tshift,3)
      call fwrite(MPI_COMM_WORLD,tshift,ng,nv,"tshift",fid_dat)
      endif
      if(is_tshift_r .eq. .true.) then
      nv = 0
      ng(1)=size(tshift_r,1)
      ng(2)=size(tshift_r,2)
      ng(3)=size(tshift_r,3)
      call fwrite(MPI_COMM_WORLD,tshift_r,ng,nv,"tshift_r",
     $fid_dat)
      endif
      if(is_u .eq. .true.) then
      nv = 0
      ng(1)=size(u,1)
      ng(2)=size(u,2)
      call fwrite(MPI_COMM_WORLD,u,ng,nv,"u",fid_dat)
      endif
      if(is_uj .eq. .true.) then
      nv = 0
      ng(1)=size(uj,1)
      ng(2)=size(uj,2)
      call fwrite(MPI_COMM_WORLD,uj,ng,nv,"uj",fid_dat)
      endif
      if(is_wgt .eq. .true.) then
      nv = 0
      ng(1)=size(wgt,1)
      call fwrite(MPI_COMM_WORLD,wgt,ng(1),nv(1),"wgt",fid_dat)
      endif
      if(is_igbs .eq. .true.) then
      nv = 0
      ng(1)=size(igbs,1)
      ng(2)=size(igbs,2)
      call fwrite(MPI_COMM_WORLD,igbs,ng,nv,"igbs",fid_dat)
      endif
      if(is_iplf_bk .eq. .true.) then
      nv = 0
      ng(1)=size(iplf_bk,1)
      ng(2)=size(iplf_bk,2)
      call fwrite(MPI_COMM_WORLD,iplf_bk,ng,nv,"iplf_bk",
     $fid_dat)
      endif
      sz_cg=size(cg)
      sz_cgr=size(cgr)
      sz_indxcg=size(indxcg)
      sz_jcg=size(jcg)
      sz_jcgr=size(jcgr)
      sz_k_list1=size(k_list,1)
      sz_k_list2=size(k_list,2)
      sz_ksym_in_npnt=size(ksym_in_npnt)

      call fwrite(MPI_COMM_WORLD,sz_cg,"sz_cg",fid_dat)
      call fwrite(MPI_COMM_WORLD,sz_cgr,"sz_cgr",fid_dat)
      call fwrite(MPI_COMM_WORLD,sz_indxcg,"sz_indxcg",fid_dat)
      call fwrite(MPI_COMM_WORLD,sz_jcg,"sz_jcg",fid_dat)
      call fwrite(MPI_COMM_WORLD,sz_jcgr,"sz_jcgr",fid_dat)
      call fwrite(MPI_COMM_WORLD,sz_k_list1,"sz_k_list1",
     $fid_dat)
      call fwrite(MPI_COMM_WORLD,sz_k_list2,"sz_k_list2",
     $fid_dat)
      call fwrite(MPI_COMM_WORLD,sz_ksym_in_npnt,
     $"sz_ksym_in_npnt",fid_dat)
      is_cg=allocated(cg)
      is_cgr=allocated(cgr)
      is_indxcg=allocated(indxcg)
      is_jcg=allocated(jcg)
      is_jcgr=allocated(jcgr)
      is_k_list=allocated(k_list)
      is_ksym_in_npnt=allocated(ksym_in_npnt)
      call fwrite(MPI_COMM_WORLD,is_cg,"is_cg",fid_dat)
      call fwrite(MPI_COMM_WORLD,is_cgr,"is_cgr",fid_dat)
      call fwrite(MPI_COMM_WORLD,is_indxcg,"is_indxcg",fid_dat)
      call fwrite(MPI_COMM_WORLD,is_jcg,"is_jcg",fid_dat)
      call fwrite(MPI_COMM_WORLD,is_jcgr,"is_jcgr",fid_dat)
      call fwrite(MPI_COMM_WORLD,is_k_list,"is_k_list",fid_dat)
      call fwrite(MPI_COMM_WORLD,is_ksym_in_npnt,
     $"is_ksym_in_npnt",fid_dat)
      if(is_cg .eq. .true.) then
      nv = 0
      ng(1)=size(cg,1)
      call fwrite(MPI_COMM_WORLD,cg,ng(1),nv(1),"cg",fid_dat)
      endif
      if(is_cgr .eq. .true.) then
      nv = 0
      ng(1)=size(cgr,1)
      call fwrite(MPI_COMM_WORLD,cgr,ng(1),nv(1),"cgr",fid_dat)
      endif
      if(is_indxcg .eq. .true.) then
      nv = 0
      ng(1)=size(indxcg,1)
      call fwrite(MPI_COMM_WORLD,indxcg,ng(1),nv(1),"indxcg",
     $fid_dat)
      endif
      if(is_jcg .eq. .true.) then
      nv = 0
      ng(1)=size(jcg,1)
      call fwrite(MPI_COMM_WORLD,jcg,ng(1),nv(1),"jcg",fid_dat)
      endif
      if(is_jcgr .eq. .true.) then
      nv = 0
      ng(1)=size(jcgr,1)
      call fwrite(MPI_COMM_WORLD,jcgr,ng(1),nv(1),"jcgr",
     $fid_dat)
      endif
      if(is_k_list .eq. .true.) then
      nv = 0
      ng(1)=size(k_list,1)
      ng(2)=size(k_list,2)
      call fwrite(MPI_COMM_WORLD,k_list,ng,nv,"k_list",fid_dat)
      endif
      if(is_ksym_in_npnt .eq. .true.) then
      nv = 0
      ng(1)=size(ksym_in_npnt,1)
      call fwrite(MPI_COMM_WORLD,ksym_in_npnt,ng(1),nv(1),
     $"ksym_in_npnt",fid_dat)
      endif
      end
