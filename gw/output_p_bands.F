      subroutine output_p_bands(key,pol)
c     key = 1 for P_0
c     key = 2 for delta_P      
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
      integer, intent(in) :: key
	complex*16, intent(in) :: pol(npb_g2_w2_max,npb_g2_w2_max,
     &	                            ndim3_nu,npnt)
	integer :: i_len,np,k,ns,ind_nu,i_nu,j,i
	real*8, allocatable :: p(:,:)
	np=min(5,npnt)
	allocate(p(0:n_nu,np))
	p=0.d0
	do k=1,np
	  ns=npb_g2_w2(k)
	  do ind_nu=1,ndim3_nu
	    i_nu=me3_tau*ndim3_nu+ind_nu-1
	    do j=1,ns
	      do i=1,ns
	        p(i_nu,k)=p(i_nu,k)+pw_p_overl_s(i,k)*pol(i,j,ind_nu,k)
     &                           *conjg(pw_p_overl_s(j,k))
            enddo
          enddo
        enddo
      enddo
	if(nproc_tau/=1) call dgop(p,np*(n_nu+1),'  +',comm_pnk)
	i_len=len_trim(allfile)
	if(maswrk) then
	  if(key==1) open(3,file=allfile(1:i_len)//'_P_0_vertex_nu')
	  if(key==2) open(3,file=allfile(1:i_len)//'_P_1_vertex_nu')
	  do i_nu=0,n_nu
 	    write(3,'(f11.5,5(1x,e12.5))')w_nu(i_nu)*evolt/2,
     &	                                (real(p(i_nu,k)),k=1,np)
	  enddo
	  close(3)
	endif
	deallocate(p)
      end