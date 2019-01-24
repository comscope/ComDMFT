      subroutine bands
****************************************************************
*   Control module for package of programs solved              *
*         eigenvalue problem of the APW method                 *
****************************************************************
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
      implicit none
	integer :: k,ind_k,ispin
	real*8 :: err_val,err_slo
	real*8, allocatable :: v_rs_mm(:,:,:)
* ****** compute eigenvalues and eigenvectors ************
	err_val=0.d0
	err_slo=0.d0
      if(ubi/='dft') then
	  allocate(v_rs_mm(n_pbmtm,n_pbmtm,ndim_nrr(me3_k+1)))
	  call v_mm_from_q(v_rs_mm)
	endif
      do ispin=1,nspin      !!! do over spins
      	if(ubi/='dft') call sigma_prelim(ispin,v_rs_mm)
        if(nspin.eq.2.and.maswrk)
     &    write(iun,*)'&&&&&&&&&   ISPIN   = ',ispin,'  &&&&&&&&&&&&'
	  n_bnd(:,ispin)=0
	  ev_bnd(:,:,:,ispin)=(0.d0,0.d0)
	  z_bnd(:,:,:,ispin)=(0.d0,0.d0)
	  e_bnd(:,:,ispin)=0.d0
	  do ind_k=1,ndim3_k(me3_k+1)
	    k=n3_mpi_k(me3_k+1)+ind_k
	    call lda_k_point(ind_k,ispin,pnt(1,k),err_val,err_slo,
     &	                   e_bnd(1,k,ispin),n_bnd(k,ispin),
     &                     z_bnd(1,1,ind_k,ispin),
     &                     ev_bnd(1,1,ind_k,ispin),indgb(1,k),nbask(k))
c     &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
          if(k.eq.1) call timel('****** LDA_K_POINT finished ********')
        enddo                !!! over k_ind
	  if(nproc_k/=1) then
	    call igop(n_bnd(1,ispin),npnt,'  +',comm_pnt)
	    call dgop(e_bnd(1,1,ispin),nbndf*npnt,'  +',comm_pnt)
	  endif
	  if(ubi/='dft') deallocate(pvv,pvg,pgg)
c -------------------------------------------------------------------
      enddo   !! over ispin
      if(ubi/='dft') deallocate(v_rs_mm)
	if(nproc_k/=1) then
	  call dgop(err_val,1,'max',comm_pnt)
	  call dgop(err_slo,1,'max',comm_pnt)
	endif
	if(maswrk) write(iun,'(a39,2e12.5)')
     &	' BANDS: Max mismatch in value, slope : ',err_val,err_slo
c      call symmetrization_z
c -------------------------------------------------------------
      do ispin=1,nspin
	  if(ubi=='dft') call output_eigen_states('_dft',0,ispin)
	  if(ubi==' hf') call output_eigen_states('_ hf',1,ispin)
      enddo   !! over ispin
      call timel('******* bands finished *************')
      end
