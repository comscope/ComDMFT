      subroutine write_gwl_bni
      use manager_mod
      use atom_mod
      use solid_mod
      use parallel_mod
      use units_mod
!      use vertex_mod
      use etot_mod
!      use models_mod
      implicit none
      include 'mpif.h'
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
        open(887, file='k_tau_freq.dat')
        write(887, '(4i10)') npnt, n_tau, n_omega,n_nu
        close(887)
       open(888,file='gwl.bni',form='unformatted')
      write(888)limlb,limlb1,limlpb,max_e_pb,maxb,maxel,
     $maxel_loc,maxlfpb,maxlfun,maxnrad,maxpb,maxpbr,maxntle,maxwf,
     $max_rl,n_pbmt,n_pbmtm,n_pbtot,nfun,nlb2,nrmax,e_c_solid_pw,
     $e_x_solid_pw,allfile,betta_t,chem_pot,chem_pot_qp,complex_ro,
     $eps_pb,interp_nu_d,irel,itest,magn,maxylm,n_omega,n_nu,n_nu_d,
     $n_nu_exa,n_tau,nrel,nrel_core,nspin,nspin_0,nspin_1,ubi,amega,
     $gbas,k_star_max,max_r_star,maxwig,maxwig_j,metal,n_symb,natom,
     $nbasmpw,nbndf,ndiv,ndtau,nkd,nkr,nplw_gw,nplw_rm,npnt,nqdiv,
     $nr_full,nr_red,nrdiv,nrr_max,nrr_max_cor,nsort,numplw,par,r0b0,
     $rb0,rbas,symmorph,v_gg_00,ngroup,maxplw,nplwdif,nbndtm,maxia,
     $maxib,maxic,nplwbas,clight
      endif
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
      if(maswrk) then
      write(888)is_augm,is_dr,is_gfun,is_gfund,is_h,is_ind_wf,
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
      if(is_augm .eq. .true.) then
       write(888)augm
       endif
      if(is_dr .eq. .true.) then
       write(888)dr
       endif
      if(is_gfun .eq. .true.) then
       write(888)gfun
       endif
      if(is_gfund .eq. .true.) then
       write(888)gfund
       endif
      if(is_h .eq. .true.) then
       write(888)h
       endif
      if(is_ind_wf .eq. .true.) then
       write(888)ind_wf
       endif
      if(is_indbasa .eq. .true.) then
       write(888)indbasa
       endif
      if(is_indfun0 .eq. .true.) then
       write(888)indfun0
       endif
      if(is_io_lem .eq. .true.) then
       write(888)io_lem
       endif
      if(is_iob2 .eq. .true.) then
       write(888)iob2
       endif
      if(is_iopb .eq. .true.) then
       write(888)iopb
       endif
      if(is_lf_isz .eq. .true.) then
       write(888)lf_isz
       endif
      if(is_lfun .eq. .true.) then
       write(888)lfun
       endif
      if(is_lfun_pb .eq. .true.) then
       write(888)lfun_pb
       endif
      if(is_lfunm .eq. .true.) then
       write(888)lfunm
       endif
      if(is_li_val .eq. .true.) then
       write(888)li_val
       endif
      if(is_lim_pb_mt .eq. .true.) then
       write(888)lim_pb_mt
       endif
      if(is_lm_isz .eq. .true.) then
       write(888)lm_isz
       endif
      if(is_lmb .eq. .true.) then
       write(888)lmb
       endif
      if(is_lmpb .eq. .true.) then
       write(888)lmpb
       endif
      if(is_n_pbmt0 .eq. .true.) then
       write(888)n_pbmt0
       endif
      if(is_nrad .eq. .true.) then
       write(888)nrad
       endif
      if(is_ntle .eq. .true.) then
       write(888)ntle
       endif
      if(is_ntle_pb .eq. .true.) then
       write(888)ntle_pb
       endif
      if(is_r .eq. .true.) then
       write(888)r
       endif
      if(is_r_atom .eq. .true.) then
       write(888)r_atom
       endif
      if(is_rl .eq. .true.) then
       write(888)rl
       endif
      if(is_smt .eq. .true.) then
       write(888)smt
       endif
      if(is_txtel .eq. .true.) then
       write(888)txtel
       endif
      if(is_bett .eq. .true.) then
       write(888)bett
       endif
      if(is_et_sum_nu .eq. .true.) then
       write(888)et_sum_nu
       endif
      if(is_ind_nu_int .eq. .true.) then
       write(888)ind_nu_int
       endif
      if(is_indxcgr .eq. .true.) then
       write(888)indxcgr
       endif
      if(is_tau_mesh .eq. .true.) then
       write(888)tau_mesh
       endif
      if(is_w_nu .eq. .true.) then
       write(888)w_nu
       endif
      if(is_w_nu_d .eq. .true.) then
       write(888)w_nu_d
       endif
      if(is_e_bnd .eq. .true.) then
       write(888)e_bnd
       endif
      if(is_e_qp .eq. .true.) then
       write(888)e_qp
       endif
      if(is_g_sym_0 .eq. .true.) then
       write(888)g_sym_0
       endif
      if(is_gbs .eq. .true.) then
       write(888)gbs
       endif
      if(is_gbs_sym .eq. .true.) then
       write(888)gbs_sym
       endif
      if(is_i_kref .eq. .true.) then
       write(888)i_kref
       endif
      if(is_ifia .eq. .true.) then
       write(888)ifia
       endif
      if(is_ig_pair .eq. .true.) then
       write(888)ig_pair
       endif
      if(is_ind_k_star_back .eq. .true.) then
       write(888)ind_k_star_back
       endif
      if(is_ind_r_star .eq. .true.) then
       write(888)ind_r_star
       endif
      if(is_ind_r_star_back .eq. .true.) then
       write(888)ind_r_star_back
       endif
      if(is_indgb .eq. .true.) then
       write(888)indgb
       endif
      if(is_indplw .eq. .true.) then
       write(888)indplw
       endif
      if(is_indpw_gw .eq. .true.) then
       write(888)indpw_gw
       endif
      if(is_indpw_r .eq. .true.) then
       write(888)indpw_r
       endif
      if(is_ip .eq. .true.) then
       write(888)ip
       endif
      if(is_ip_k .eq. .true.) then
       write(888)ip_k
       endif
      if(is_ip_r .eq. .true.) then
       write(888)ip_r
       endif
      if(is_iplf_gk .eq. .true.) then
       write(888)iplf_gk
       endif
      if(is_ir_ref .eq. .true.) then
       write(888)ir_ref
       endif
      if(is_is .eq. .true.) then
       write(888)is
       endif
      if(is_itshift_r .eq. .true.) then
       write(888)itshift_r
       endif
      if(is_k_group .eq. .true.) then
       write(888)k_group
       endif
      if(is_k_star .eq. .true.) then
       write(888)k_star
       endif
      if(is_k_sym_0 .eq. .true.) then
       write(888)k_sym_0
       endif
      if(is_list_sym_k .eq. .true.) then
       write(888)list_sym_k
       endif
      if(is_n_bnd .eq. .true.) then
       write(888)n_bnd
       endif
      if(is_nbask .eq. .true.) then
       write(888)nbask
       endif
      if(is_nhsort .eq. .true.) then
       write(888)nhsort
       endif
      if(is_nplw_r .eq. .true.) then
       write(888)nplw_r
       endif
      if(is_nplwgw .eq. .true.) then
       write(888)nplwgw
       endif
      if(is_nrr_red .eq. .true.) then
       write(888)nrr_red
       endif
      if(is_num_r_star .eq. .true.) then
       write(888)num_r_star
       endif
      if(is_num_sym_k .eq. .true.) then
       write(888)num_sym_k
       endif
      if(is_pnt .eq. .true.) then
       write(888)pnt
       endif
      if(is_r0_pair .eq. .true.) then
       write(888)r0_pair
       endif
      if(is_r_group .eq. .true.) then
       write(888)r_group
       endif
      if(is_r_pnt .eq. .true.) then
       write(888)r_pnt
       endif
      if(is_r_pnt_ind .eq. .true.) then
       write(888)r_pnt_ind
       endif
      if(is_shift .eq. .true.) then
       write(888)shift
       endif
      if(is_sovi .eq. .true.) then
       write(888)sovi
       endif
      if(is_sovr .eq. .true.) then
       write(888)sovr
       endif
      if(is_sym_label .eq. .true.) then
       write(888)sym_label
       endif
      if(is_tau .eq. .true.) then
       write(888)tau
       endif
      if(is_tshift .eq. .true.) then
       write(888)tshift
       endif
      if(is_tshift_r .eq. .true.) then
       write(888)tshift_r
       endif
      if(is_u .eq. .true.) then
       write(888)u
       endif
      if(is_uj .eq. .true.) then
       write(888)uj
       endif
      if(is_wgt .eq. .true.) then
       write(888)wgt
       endif
      if(is_igbs .eq. .true.) then
       write(888)igbs
       endif
      if(is_iplf_bk .eq. .true.) then
       write(888)iplf_bk
       endif
      endif
      sz_cg=size(cg)
      sz_cgr=size(cgr)
      sz_indxcg=size(indxcg)
      sz_jcg=size(jcg)
      sz_jcgr=size(jcgr)
      sz_k_list1=size(k_list,1)
      sz_k_list2=size(k_list,2)
      sz_ksym_in_npnt=size(ksym_in_npnt)

      if(maswrk) then
      write(888)sz_cg,sz_cgr,sz_indxcg,sz_jcg,sz_jcgr,
     $sz_k_list1,sz_k_list2,sz_ksym_in_npnt
      endif
      is_cg=allocated(cg)
      is_cgr=allocated(cgr)
      is_indxcg=allocated(indxcg)
      is_jcg=allocated(jcg)
      is_jcgr=allocated(jcgr)
      is_k_list=allocated(k_list)
      is_ksym_in_npnt=allocated(ksym_in_npnt)
      if(maswrk) then
      write(888)is_cg,is_cgr,is_indxcg,is_jcg,is_jcgr,is_k_list,
     $is_ksym_in_npnt
      if(is_cg .eq. .true.) then
       write(888)cg
       endif
      if(is_cgr .eq. .true.) then
       write(888)cgr
       endif
      if(is_indxcg .eq. .true.) then
       write(888)indxcg
       endif
      if(is_jcg .eq. .true.) then
       write(888)jcg
       endif
      if(is_jcgr .eq. .true.) then
       write(888)jcgr
       endif
      if(is_k_list .eq. .true.) then
       write(888)k_list
       endif
      if(is_ksym_in_npnt .eq. .true.) then
       write(888)ksym_in_npnt
       endif
      endif
      is_ab_nu_from_tau=allocated(ab_nu_from_tau)
      is_tau_from_nu=allocated(tau_from_nu)
      is_ev_bnd=allocated(ev_bnd)
      is_q_qp=allocated(q_qp)
      is_s_it_r0=allocated(s_it_r0)
      is_v_ii=allocated(v_ii)
      is_v_mm=allocated(v_mm)
      is_w_rs_mi=allocated(w_rs_mi)
      is_z_bnd=allocated(z_bnd)
      is_ev_qp=allocated(ev_qp)
      is_z_qp=allocated(z_qp)
      if(maswrk) then
      write(888)is_ab_nu_from_tau,is_tau_from_nu,is_ev_bnd,
     $is_q_qp,is_s_it_r0,is_v_ii,is_v_mm,is_w_rs_mi,is_z_bnd,is_ev_qp,
     $is_z_qp
      close(888)
      endif
      i_len=len_trim(allfile)

       if (me3_tau .eq. 0) then
        do ind_k=1,ndim3_kk(me3_k+1)
          k0=n3_mpi_kk(me3_k+1)+ind_k
          if (k0 .lt. 10) then
            write(k0_ch, '(i1)') k0
          elseif (k0 .lt. 100) then
            write(k0_ch, '(i2)') k0
          elseif (k0 .lt. 1000) then
            write(k0_ch, '(i3)') k0
          elseif (k0 .lt. 10000) then
            write(k0_ch, '(i4)') k0
          endif
      if(is_s_it_r0 .eq. .true.) then
      open(3,file=allfile(1:i_len)//'_s_it_r0_'//trim(k0_ch)
     $//'_'//trim(ubi)//'.uca',form='unformatted')
      write(3)s_it_r0(:,:,ind_k)
      close(3)
      endif
      if(is_w_rs_mi .eq. .true.) then
      open(3,file=allfile(1:i_len)//'_w_rs_mi_'//trim(k0_ch)
     $//'_'//trim(ubi)//'.uca',form='unformatted')
      write(3)w_rs_mi(:,:,ind_k)
      close(3)
      endif
      enddo
      endif

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
      open(3,file=allfile(1:i_len)//'_ev_bnd_'//trim(k0_ch)
     $//'_'//trim(ubi)//'.uca',form='unformatted')
      do ispin=1, nspin_0
      write(3)ev_bnd(:,:,ind_k,ispin)
      enddo
      close(3)
      endif
      if(is_q_qp .eq. .true.) then
      open(3,file=allfile(1:i_len)//'_q_qp_'//trim(k0_ch)
     $//'_'//trim(ubi)//'.uca',form='unformatted')
      do ispin=1, nspin_0
      write(3)q_qp(:,:,ind_k,ispin)
      enddo
      close(3)
      endif
      if(is_v_ii .eq. .true.) then
      open(3,file=allfile(1:i_len)//'_v_ii_'//trim(k0_ch)
     $//'_'//trim(ubi)//'.uca',form='unformatted')
      write(3)v_ii(:,:,ind_k)
      close(3)
      endif
      if(is_v_mm .eq. .true.) then
      open(3,file=allfile(1:i_len)//'_v_mm_'//trim(k0_ch)
     $//'_'//trim(ubi)//'.uca',form='unformatted')
      write(3)v_mm(:,:,ind_k)
      close(3)
      endif
      if(is_z_bnd .eq. .true.) then
      open(3,file=allfile(1:i_len)//'_z_bnd_'//trim(k0_ch)
     $//'_'//trim(ubi)//'.uca',form='unformatted')
      do ispin=1, nspin_0
      write(3)z_bnd(:,:,ind_k,ispin)
      enddo
      close(3)
      endif
      if(is_ev_qp .eq. .true.) then
      open(3,file=allfile(1:i_len)//'_ev_qp_'//trim(k0_ch)
     $//'_'//trim(ubi)//'.uca',form='unformatted')
      do ispin=1, nspin_0
      write(3)ev_qp(:,:,ind_k,ispin)
      enddo
      close(3)
      endif
      if(is_z_qp .eq. .true.) then
      open(3,file=allfile(1:i_len)//'_z_qp_'//trim(k0_ch)
     $//'_'//trim(ubi)//'.uca',form='unformatted')
      do ispin=1, nspin_0
      write(3)z_qp(:,:,ind_k,ispin)
      enddo
      close(3)
      endif
        enddo
       endif
      i_len=len_trim(allfile)
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
      open(3,file=allfile(1:i_len)//'_ab_nu_from_tau_'//trim(tau_ch)
     $//'_'//trim(ubi)//'.uca',form='unformatted')
      write(3)ab_nu_from_tau(ind_tau,:,:)
      close(3)
      endif
      if(is_tau_from_nu .eq. .true.) then
      open(3,file=allfile(1:i_len)//'_tau_from_nu_'//trim(tau_ch)
     $//'_'//trim(ubi)//'.uca',form='unformatted')
      write(3)tau_from_nu(:,ind_tau,:)
      close(3)
      endif
        enddo
       endif
      i_len=len_trim(allfile)
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
      i_len=len_trim(allfile)
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
      i_len=len_trim(allfile)
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

      end
