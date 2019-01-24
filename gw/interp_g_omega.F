      subroutine interp_g_omega(ab_in,ab_cheb,f_out,x_out,n_out,m,n,key)
c ------- We assume that G_ij(-w) = G_ij*(w), i.e. G(t) - real -----
c     key = 0 - for Gx+Gc
c     key = 1 - only for Gc
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
      integer, intent(in) :: n_out,n,m,key
	real*8, intent(in) :: ab_in(m,n,2,0:n_omega),x_out(0:n_out),
     &                      ab_cheb(m,n,2,0:n_omega_asy-2)
	complex*16, intent(out) :: f_out(m,n,0:n_out)
	integer :: i,j,io
	real*8 :: w,w1
	real*8, allocatable :: aa(:,:,:)
	allocate(aa(m,n,2))
	do io=0,n_out
	  w1=x_out(io)
	  w=abs(w1)
	  call ferm_interp_omega(m*n,ab_in,ab_cheb,w,aa,key)
	  do j=1,n
	    do i=1,m
	      f_out(i,j,io)=dcmplx(aa(i,j,1),aa(i,j,2))
	    enddo
	  enddo
	  if(w1<0) f_out(:,:,io)=conjg(f_out(:,:,io))
	enddo   !! over io
	deallocate(aa)
      end