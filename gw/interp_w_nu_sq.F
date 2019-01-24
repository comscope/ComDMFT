      subroutine interp_w_nu_sq(f_in,chb,f_out,x_out,n_out,n)
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
      integer, intent(in) :: n_out,n
	real*8, intent(in) :: f_in(n,n,0:n_nu),chb(n,n,0:n_nu_asy-2),
     &	                    x_out(0:n_out)
	complex*16, intent(out) :: f_out(n,n,0:n_out)
	integer :: io,nn
	real*8, allocatable :: aa(:,:)
	nn=n*n
	allocate(aa(n,n))
	do io=0,n_out
	  call boson_interp_nu_v(nn,f_in,chb,x_out(io),aa)
	  call unpack_hermit(f_out(1,1,io),aa,n,n,n,(0.d0,0.d0),
     &	                   (1.d0,0.d0))
	enddo   !! over io
	deallocate(aa)
      end