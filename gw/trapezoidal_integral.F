      real*8 function trapezoidal_integral(y,x,n)
      implicit none
	integer, intent(in) :: n
      real*8, intent(in) :: y(0:n),x(0:n)
	integer :: i
	real*8 :: s
	s=0.d0
	do i=1,n
	  s=s+0.5d0*(y(i-1)+y(i))*(x(i)-x(i-1))
	enddo
	trapezoidal_integral=s
      end