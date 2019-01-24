      subroutine etot_gw_2
c -------- The OUTPUT of GW ENERGY -------------------------------------
	use atom_mod
	use etot_mod
	use manager_mod
	use parallel_mod
	use solid_mod
      implicit none
	real*8 :: gvh,t_kin_x,t_kin_c,e_coul,t_kin,e_x,e_c,e_xc,etot,
     &	        e_xc_pw,etot_pw,e_x_pw,ts,psi_pw,mu_n,g_rpa_lda,
     &          mu_n_dft,e_rpa_lda
      gvh=0.0d0
      if(nproc_k/=1) then
	  call DGOP(gx_vh_val,1,'  +',comm_pnt)
	  call DGOP(ebnd_x,1,'  +',comm_pnt)
	  call DGOP(gx_sx_val,1,'  +',comm_pnt)
!         call DGOP(gx_vxc_val,1,'  +',comm_pnt) ! redundant
	  call DGOP(e_magn,1,'  +',comm_pnt)
	  call DGOP(d_free_x,1,'  +',comm_pnt)
	  if(ubi==' qp'.or.ubi==' gw'.or.ubi=='gwg') then
	    call DGOP(gc_vh_val,1,'  +',comm_pnt)
	    call DGOP(ebnd_c,1,'  +',comm_pnt)
	    call DGOP(gc_sx_val,1,'  +',comm_pnt)
	    call DGOP(gc_sc_val,1,'  +',comm_pnt)
	    call DGOP(gx_sc_val,1,'  +',comm_pnt)
	    call DGOP(gc_sc_val_w,1,'  +',comm_pnt)
	    call DGOP(gx_sc_val_w,1,'  +',comm_pnt)
	    call DGOP(gc_vxc_val,1,'  +',comm_pnt)
	    call DGOP(d_free_c,1,'  +',comm_pnt)
	  endif
	endif
      if(nproc_tau/=1) then
	  if(ubi_0=='dft') 
     &	  call brdcst(msgdbl,gc_vxc_val,8,0,comm_pnk)
	  if(ubi==' qp'.or.ubi==' gw'.or.ubi=='gwg') then
	    call brdcst(msgdbl,gc_vh_val,8,0,comm_pnk)
	    call brdcst(msgdbl,gc_sx_val,8,0,comm_pnk)
	    call brdcst(msgdbl,ebnd_c,8,0,comm_pnk)
	    call DGOP(gc_sc_val,1,'  +',comm_pnk)
	    call DGOP(gx_sc_val,1,'  +',comm_pnk)
	    call DGOP(gc_sc_val_w,1,'  +',comm_pnk)
	    call DGOP(gx_sc_val_w,1,'  +',comm_pnk)
	  endif
	endif
	if(maswrk)
     &write(iun,139) z_vnucl,ro_vh_old,ro_vh_new,ecor,ebnd_x,ebnd_c,
     &               gx_vxc_cor/2,gx_vxc_val/2,gc_vxc_val/2,
     &               gx_vh_cor/2,gx_vh_val/2,gc_vh_val/2,
     &               gx_sx_cor/2,gx_sx_val/2,gc_sx_val/2,
     &               gx_sc_val/2,gc_sc_val/2,gx_sc_val_w/2,
     &               gc_sc_val_w/2
  139 format(/' Z_VNUCL=',f17.8,' ROVHOLD=',f17.8,' ROVHNEW=',f17.8,
     &      /' ECOR=   ',f16.8,'   EBND_X= ',f15.8,'  EBND_C = ',f15.8,
     &     /'GX_VXC_COR=',f15.8,' GX_VXC_VAL=',f14.8,' GC_VXC_V=',f14.8,
     &      /'GX_VH_COR=',f16.8,' GX_VH_VAL=',f15.8,' GC_VH_VAL=',f15.8,
     &      /'GX_SX_COR=',f15.8,' GX_SX_VAL=',f15.8,' GC_SX_VAL=',f15.8,
     &      /' GX_SC_VAL = ',f15.8,' GC_SC_VAL = ',f15.8,
     &      /' GX_SC_VAL_W=',f15.8,' GC_SC_VAL_W=',f15.8)
      mu_n=chem_pot*(sum(zcor)+nelec)
      mu_n_dft=chem_pot_lda*(sum(zcor)+nelec)
      t_kin_c=0.d0
	if(ubi=='dft') then
	  e_coul=z_vnucl+0.5d0*ro_vh_new
	else
	  gvh=gx_vh_cor+gx_vh_val+gc_vh_val
	  e_coul=z_vnucl+0.5d0*gvh
	endif
	t_kin_x=ecor+ebnd_x
	if(ubi=='dft') then
	  t_kin_x=t_kin_x-ro_vh_old-vxch
	else 
	  t_kin_x=t_kin_x-(gx_vh_cor+gx_vh_val)
	  t_kin_x=t_kin_x-gx_sx_cor-gx_sx_val
	endif
	if(ubi==' qp'.or.ubi==' gw'.or.ubi=='gwg') then
	  t_kin_c=ebnd_c-gc_vh_val-gc_sx_val
	endif
c -------------- Convolutions SIGMA*GREEN -----------------------
	e_x=0.d0
	e_c=0.d0
	d_free_x=d_free_x+ecor-chem_pot*sum(zcor)
	if(ubi=='dft') then
	  e_xc=exch_dft
	  grand_pot=d_free_x-ro_vh_old-vxch
	  grand_pot=grand_pot+e_coul+e_xc+e_b_ext
	  grand_pot_dft=grand_pot
	else
	  e_x=0.5d0*(gx_sx_cor+gx_sx_val+gc_sx_val)
	  grand_pot_x=d_free_x+e_coul-gvh-e_x
	  grand_pot=grand_pot_x
	  if(ubi==' qp'.or.ubi==' gw'.or.ubi=='gwg') then
	    e_c=0.5d0*(gx_sc_val+gc_sc_val)
	    if(ubi=='gwg') then
            psi_pw=e_vertex_p_w1-e_vertex_pw_tau0
	      phi_functional_c=phi_functional_c+psi_functional_vrt
     &	                    -psi_pw
	    endif
	    grand_pot=grand_pot_x+phi_functional_c-2.d0*e_c+d_free_c
	  endif
	  e_xc=e_x+e_c
	endif
	t_kin=t_kin_x+t_kin_c
	etot=t_kin+e_coul+e_xc+e_b_ext
	if(ubi=='dft') then
	  e_tot_lda=etot
	  exc_tot_lda=e_xc
	  ts=ebnd_x-d_free_x-mu_n
	else
	  ts=etot-grand_pot-mu_n
	endif
	if(iter==iter_max_scf.and.iter==iter_dft+1) then
c ---------- RPA(LDA) approach --------------------------------
        g_rpa_lda=grand_pot_dft+e_x_dft+phi_functional_dft-exch_dft
        e_rpa_lda=e_tot_lda+e_x_dft+phi_functional_dft-exch_dft
	  if(maswrk) then
	    write(iun,*)' Phi_c_RPA_DFT ',phi_functional_dft
	    write(iun,*)' Grand Potential RPA(DFT)',g_rpa_lda
	    write(iun,*)' Free Energy RPA(DFT) ',g_rpa_lda+mu_n_dft
	    write(iun,*)' Internal Energy RPA(DFT) ',e_rpa_lda
	  endif
	endif
	if(maswrk) then
	  write(iun,*)' Energy of valence Spin-Orbit Interaction ',
     &	                   e_spin_orbit
	  write(iun,*)' Energy of full Spin-Orbit Interaction ',
     &	                   e_spin_orbit+e_spin_orbit_core
	  write(iun,*)' D_Free_X, D_Free_C ',d_free_x,d_free_c
	  write(iun,*)' Phi_c_RPA ',phi_functional_c
	  write(iun,*)' Grand Potential ',grand_pot
	  write(iun,*)' Free Energy F ',grand_pot+mu_n
	  write(iun,*)' Entropy (Ry/K) ',ts/temperature
	  if(ubi=='gwg') then
	    write(iun,*)' E_vertex_sig_g_FULL_BND ',e_vertex_sig_g_bnd
	    write(iun,*)' E_vertex_P_W_FULL_PB    ',e_vertex_p_w_pb
	    write(iun,*)' Psi_c_VRT_Sig_G ',psi_functional_vrt
	    write(iun,*)' Psi_c_VRT_P_W ',psi_pw
	  endif
	  if(ubi==' gw'.or.ubi==' qp') then
	    e_x_pw=0.5d0*gx_sx_cor+e_x_solid_pw
	    e_xc_pw=e_x_pw+e_c_solid_pw
	    etot_pw=t_kin+e_coul+e_x+e_c_solid_pw+e_b_ext
          write(iun,140) e_coul,t_kin,e_x,e_c,e_xc,e_x_solid_pw,
     &                   e_c_solid_pw,e_xc_pw,etot,etot_pw
	    print '(a3,1x,a20,3(1x,f17.8))',ubi,': E(SG), E(PW), F = ',
     &	                             etot,etot_pw,grand_pot+mu_n
	  else
          write(iun,141) e_coul,t_kin,e_x,e_c,e_xc,etot
	    print '(a4,a11,2f22.12)',ubi,': ETOT F = ',
     &                               etot,grand_pot+mu_n
	  endif
	endif
      call timel('***** ETOT_GW_2 finished ***********')
  140 format(/' E_COUL  = ',f17.8,' T_KIN  = ',f15.8
     .       /'GS : E_X = ',f15.8,'   E_C  = ',f15.8,'  E_XC   = ',f15.8
     .       /'PW : E_X = ',f15.8,'   E_C  = ',f15.8,'  E_XC   = ',f15.8
     .       /' ETOT (GS)  = ',f17.8
     .       /' ETOT (PW)  = ',f17.8)
  141 format(/' E_COUL  = ',f17.8,' T_KIN  = ',f15.8
     .       /'GS : E_X = ',f15.8,'   E_C  = ',f15.8,'  E_XC   = ',f15.8
     .       /' ETOT (GS)  = ',f17.8)
      end
