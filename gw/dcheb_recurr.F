	subroutine dcheb_recurr(x,cheb,n)
******************************************************************
*   Calculate the Chebyshev polinomials up to order n            *
******************************************************************
      implicit none
	integer, intent(in) :: n
	real*8, intent(in) :: x
	real*8, intent(out) :: cheb(0:n)
	integer :: i
	real*8 :: x2
	x2=x+x
	cheb(0)=1.d0
	if(n.eq.0) return
	cheb(1)=x
	if(n.eq.1) return
	do i=1,n-1
	  cheb(i+1)=x2*cheb(i)-cheb(i-1)
	enddo
      end