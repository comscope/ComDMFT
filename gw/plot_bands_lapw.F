      subroutine plot_bands_lapw
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
	implicit none
	integer :: ispin,k,ind_k,nbn
	real*8 :: err_val,err_slo
	complex*16, allocatable :: zbn(:,:),evbn(:,:)
c ---------------------------------------------------------------------
	allocate(e_plot(nbndf,num_k_plot,nspin))
      allocate(zbn(nfun,nbndf))
      allocate(evbn(nbasmpw,nbndf))
	e_plot=0.d0
	nbnd_plot=100000
	err_val=0.d0
	err_slo=0.d0
	pv=0.d0
	if(irel==2) pvj=0.d0
	call vslli(2,2,v_mt_h,v_mt_xc,v_intr_h,v_intr_xc)
     &	 ! Hartree + XC
	do ispin=1,nspin
	  do ind_k=1,ndim_kplot(me3_k+1)
	    k=n_mpi_kplot(me3_k+1)+ind_k
	    call lda_k_point(k,ispin,k_plot(1,k),err_val,err_slo,
     &	                   e_plot(1,k,ispin),nbn,zbn,evbn,
     &                     indgb_plot(1,k),nbask_plot(k))
	    nbnd_plot=min(nbnd_plot,nbn)
	  enddo  !! over ind_k
	enddo    !! over ispin
	k=nbndf*num_k_plot*nspin
	if(nproc_k/=1) then
	  call DGOP(e_plot,k,'  +',comm_pnt)
	  call IGOP(nbnd_plot,1,'min',comm_pnt)
	endif
	call bands_gnuplot(nbnd_plot,nbndf,e_plot)
	deallocate(e_plot,zbn,evbn)
      call timel('***** PLOT_BANDS_LAPW finished *****')
      end