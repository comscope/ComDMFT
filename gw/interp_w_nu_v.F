      subroutine interp_w_nu_v(f_in,chb,f_out,x_out,n_out,n)
c ------- We assume that W_ij(-w) = W_ij(w),  -----
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
      integer, intent(in) :: n_out,n
	real*8, intent(in) :: f_in(n*(n+1)/2,0:n_nu),
     &	                    chb(n*(n+1)/2,0:n_nu_asy-2),x_out(0:n_out)
	real*8, intent(out) :: f_out(n,n,0:n_out)
	integer :: i,j,io,nn,ij
	real*8, allocatable :: aa(:)
	nn=n*(n+1)/2
	allocate(aa(nn))
	do io=0,n_out
	  call boson_interp_nu_v(nn,f_in,chb,x_out(io),aa)
	  ij=0
	  do j=1,n
	    do i=1,j
	      ij=ij+1
	      f_out(i,j,io)=aa(ij)
	      f_out(j,i,io)=aa(ij)
	    enddo
	  enddo
	enddo   !! over io
	deallocate(aa)
      end