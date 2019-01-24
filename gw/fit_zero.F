      subroutine fit_zero(f,r)
      implicit none
	real*8, intent(inout) :: f(0:3),r(0:3)
	integer :: i
	real*8 :: a(3),b(3),c(3),det3,det_all
	a=1.d0
	b=r(1:)
	do i=1,3
	  c(i)=r(i)**2
	enddo
	det_all=det3(a,b,c)
	f(0)=det3(f(1),b,c)/det_all
      end