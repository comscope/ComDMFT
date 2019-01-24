      real*8 function conv_hh(a,b,n,n0,n1)
c	calculates Sum_ij  a_ij*b_ji for two hermition matrices
c     They are stored in R*8 square form
      implicit none
	integer, intent(in) :: n,n0,n1
	real*8, intent(in) :: a(n0,n0),b(n1,n1)
      integer :: i
	conv_hh=0.0
	do i=1,n
	  conv_hh=conv_hh+2.d0*(dot_product(a(1:i-1,i),b(1:i-1,i))
     &	                     +dot_product(a(i,1:i-1),b(i,1:i-1)))
	  conv_hh=conv_hh+a(i,i)*b(i,i)
	enddo
      end