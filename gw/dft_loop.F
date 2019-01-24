      subroutine dft_loop
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	implicit none
	integer :: ispin,ind_k,k
	if(ubi/='dft') then
	  v_mt_h=adm_gw*v_mt_h+(1.d0-adm_gw)*v_mt_h_old
	endif
	v_mt_h_old=v_mt_h
      call rad_eqs
	pv=0.d0
	if(irel==2) pvj=0.d0
c ----------- Full Potential; L=0 term is absent ---------------
      if(ubi=='dft') then
	  call vslli(2,2,v_mt_h,v_mt_xc,v_intr_h,v_intr_xc)
      else
	  call vslli(0,2,v_mt_h,v_mt_xc,v_intr_h,v_intr_xc)
	endif
      call timel('******* vslli finished *************')
	if(iter>iter_dft.and.we_need_v) then
	  call prod_bas
	  call trans_pb_pw
	endif
      call bands
      if(iter==1) chem_pot=minval(e_bnd(n_semic+1,:,1:nspin))+0.5d0
c -----------------------------------------------------------------
	call search_mu_0(chem_pot,nelec,chem_pot)
      call timel('**** Search_mu finished ************')
c -----------------------------------------------------------------
	call set_g
      call timel('******** Set_G finished ************')
c -----------------------------------------------------------------
      call seteny
c -----------------------------------------------------------------      
      if(iter>=iter_dft.and.iter_dft<iter_max_scf) then
c ------ Getting the matrix of Hartree potential  ----
	  pv=0.d0
	  if(irel==2) pvj=0.d0
	  call vslli(0,1,v_mt_h,v_mt_xc,v_intr_h,v_intr_xc) ! only Hartree
	  do ispin=1,nspin
	    do ind_k=1,ndim3_k(me3_k+1)
	      k=n3_mpi_k(me3_k+1)+ind_k
	      call v_bloch(v_h_bb(1,1,ind_k,ispin),pnt(1,k),ispin,
     &                   n_bnd(k,ispin),nbndf,nbask(k),indgb(1,k),
     &                   ev_bnd(1,1,ind_k,ispin),z_bnd(1,1,ind_k,ispin),
     &                   pv,pvj,war,wari)
	    enddo
	  enddo
	endif
      if(iter_dft<iter_max_scf.and.iter==iter_dft) then
	  pv=0.d0
	  if(irel==2) pvj=0.d0
        call vslli(1,1,v_mt_h,v_mt_xc,v_intr_h,v_intr_xc) ! only XC
	  do ispin=1,nspin
	    do ind_k=1,ndim3_k(me3_k+1)
	      k=n3_mpi_k(me3_k+1)+ind_k
	      call v_bloch(v_xc_bb(1,1,ind_k,ispin),pnt(1,k),ispin,
     &                   n_bnd(k,ispin),nbndf,nbask(k),indgb(1,k),
     &                   ev_bnd(1,1,ind_k,ispin),z_bnd(1,1,ind_k,ispin),
     &                   pv,pvj,war,wari)
	    enddo  
	  enddo
      endif      
      if(iter_dft<iter_max_scf.and.iter>=iter_dft) then
        v_h_bb_old=v_h_bb
	  if(ubi=='dft') then
	    v_xc_bb_old=v_xc_bb
	  else
	    sigx_solid_k_old=sigx_solid_k
	  endif
	endif
      end
