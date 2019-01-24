	subroutine d_cheb_dx(x,dcheb,cheb,n)
******************************************************************
*   Calculate the derivatives of Chebyshev polinomials with respect to x
******************************************************************
      implicit none
	integer, intent(in) :: n
	real*8, intent(in) :: x,cheb(0:n)
	real*8, intent(out) :: dcheb(0:n)
	integer :: i
	real*8 :: xx
	if(x==-1.d0) then
	  do i=0,n
	    dcheb(i)=dfloat(-(-1)**i*i*i)
	  enddo
	else if(x==1.d0) then
	  do i=0,n
	    dcheb(i)=dfloat(i*i)
	  enddo
	else
	  xx=1.d0-x*x
	  dcheb(0)=0.d0
	  do i=1,n
	    dcheb(i)=dfloat(i)*(cheb(i-1)-x*cheb(i))/xx
	  enddo
	endif
      end