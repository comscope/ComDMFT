      subroutine output_vh_sigx_band
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer :: k,i_len,n,ispin,ind_k,i,ib,n0,nmin
	real*8, allocatable :: z_out(:,:,:)
	n0=nbndf
	nmin=minval(n_bnd)
	n=min(npnt,6)
	allocate(z_out(nmin,n,2))
	if(maswrk) then
	  i_len=len_trim(allfile)
	  if(ubi==' gw') then
	    open(3,file=allfile(1:i_len)//'_V_Hartree_band_k.gw')
	    open(4,file=allfile(1:i_len)//'_Sigma_X_band_k.gw')
	  else if(ubi==' hf') then
	    open(3,file=allfile(1:i_len)//'_V_Hartree_band_k.hf')
	    open(4,file=allfile(1:i_len)//'_Sigma_X_band_k.hf')
	  else if(ubi==' qp') then
	    open(3,file=allfile(1:i_len)//'_V_Hartree_band_k.qp')
	    open(4,file=allfile(1:i_len)//'_Sigma_X_band_k.qp')
	  else if(ubi=='gwg') then
	    open(3,file=allfile(1:i_len)//'_V_Hartree_band_k.gwg')
	    open(4,file=allfile(1:i_len)//'_Sigma_X_band_k.gwg')
	  endif
	endif
	do ispin=1,nspin
	  if(maswrk) write(3,*)' Ispin = ',ispin
	  if(maswrk) write(4,*)' Ispin = ',ispin
	  z_out=0.d0
	  do ind_k=1,ndim3_k(me3_k+1)
	    k=n3_mpi_k(me3_k+1)+ind_k
	    if(k<=n) then
	      do i=1,nmin
	        z_out(i,k,1)=v_h_bb(i,i,ind_k,ispin)
	        z_out(i,k,2)=sigx_solid_k(i,i,ind_k,ispin)
	      enddo
	    endif
	  enddo
c --------- Collect from all processes  ------------------------------------------
	  if(nproc_k/=1) call dgop(z_out,nmin*n*2,'  +',comm_pnt)
	  if(maswrk) then
	    write(3,*)' K-points from ',1,' to ',n
	    write(4,*)' K-points from ',1,' to ',n
	    do ib=1,nmin
	      write(3,'(i5,6(1x,f11.6))')ib,(z_out(ib,k,1),k=1,n)
	      write(4,'(i5,6(1x,f11.6))')ib,(z_out(ib,k,2),k=1,n)
	    enddo
	  endif
	enddo
	if(maswrk) then
	  close(3)
	  close(4)
	endif
	deallocate(z_out)
      end
