      subroutine output_z_factor
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer :: k,i_len,m,n,ispin,l1,l2,m0,ind_k,i,nb,info,ib,n0
	integer, allocatable :: ipiv(:)
	real*8, allocatable :: z_out(:,:)
	complex*16, allocatable :: b(:,:),work(:,:)
	n0=nbndf
	m=5
	n=npnt/m
	if(n*m<npnt) n=n+1
	allocate(z_out(n0,npnt),b(n0,n0),ipiv(n0),work(n0,n0))
	if(maswrk) then
	  i_len=len_trim(allfile)
	  if(ubi==' gw') then
	    open(3,file=allfile(1:i_len)//'_Z_factor_band_k_x.gw')
	  else if(ubi==' qp') then
	    open(3,file=allfile(1:i_len)//'_Z_factor_band_k_x.qp')
	  else if(ubi=='gwg') then
	    open(3,file=allfile(1:i_len)//'_Z_factor_band_k_x.gwg')
	  endif
	endif
	do ispin=1,nspin
	  if(maswrk) write(3,*)' Ispin = ',ispin
	  z_out=0.d0
	  do ind_k=1,ndim3_k(me3_k+1)
	    k=n3_mpi_k(me3_k+1)+ind_k
	    nb=n_bnd(k,ispin)
	    call unpack_hermit(b,z_ren(1,1,ind_k,ispin),nb,n0,n0,
     &		                 (0.d0,0.d0),(1.d0,0.d0))
          CALL zgetrf(nb,nb,b,n0,ipiv,INFO)
          call zgetri(nb,b,n0,ipiv,work,n0**2,info)
	    do i=1,nb
	      z_out(i,k)=real(b(i,i))
	    enddo
	  enddo
c --------- Collect from all processes  ------------------------------------------
	  if(nproc_k/=1) then
	    call dgop(z_out,n0*npnt,'  +',comm_pnt)
	  endif
	  do i=1,n
          l1=m*(i-1)+1
          l2=min(m*i,npnt)
	    if(maswrk) write(3,*)' K-points from ',l1,' to ',l2
	    m0=l2-l1+1
	    if(maswrk) then
	      do ib=1,n0
	        write(3,'(i5,5(1x,f12.7))')ib,(z_out(ib,k),k=1,m0)
	      enddo
	    endif
	  enddo
	enddo
	if(maswrk) then
	  close(3)
	endif
	deallocate(z_out,b,ipiv,work)
      end
