      subroutine interp_g_omega_c(ab_in,ab_cheb,f_out,x_out,n_out,n,key)
c ------- We assume that G_ij(-w) = G_ji*(w), i.e. G(t) - Hermitian ---
c --- Input and output arrays are in compact R*8 form -----------------
c     key = 0 - for Gx+Gc
c     key = 1 - only for Gc
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
      integer, intent(in) :: n_out,n,key
	real*8, intent(in) :: ab_in(n,n,2,0:n_omega),x_out(0:n_out),
     &                      ab_cheb(n,n,2,0:n_omega_asy-2)
	complex*16, intent(out) :: f_out(n,n,0:n_out)
	integer :: io
	real*8 :: w,w1
	real*8, allocatable :: aa(:,:,:)
	allocate(aa(n,n,2))
	do io=0,n_out
	  w1=x_out(io)
	  w=abs(w1)
	  call ferm_interp_omega(n*n,ab_in,ab_cheb,w,aa,key)
	  call ferm_unpack_omega(f_out(1,1,io),aa,n,n,n)
	  if(w1<0) f_out(:,:,io)=conjg(transpose(f_out(:,:,io)))
	enddo   !! over io
	deallocate(aa)
      end