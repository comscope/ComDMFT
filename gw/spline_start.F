      subroutine spline_start(h,y,a,b,c,d,n)
      implicit none
	integer, intent(in) :: n
	real*8, intent(in) :: h(n),y(0:n)
	real*8, intent(out) :: a(n),b(n),c(n),d(n)
	real*8, allocatable :: x(:)
	integer :: i,k
	do i=2,n
	  k=i-1
	  a(k)=h(k)
	  b(k)=-2.d0*(h(k)+h(i))
	  c(k)=h(i)
	  d(k)=3.d0*((y(i)-y(k))/h(i)-(y(k)-y(k-1))/h(k))
	enddo
	allocate(x(n))
	call spline(a,b,c,d,n-1,x(2))
	c(1)=0.d0
	c(2:)=x(2:)
	do i=1,n-1
	  b(i)=(y(i)-y(i-1))/h(i)-h(i)*(c(i+1)+2.d0*c(i))/3.d0
	  d(i)=(c(i+1)-c(i))/3.d0/h(i)
	  a(i)=y(i-1)
	enddo
	b(n)=(y(n)-y(n-1))/h(i)-2.d0*h(n)*c(n)/3.d0
	d(n)=-c(n)/3.d0/h(n)
	a(n)=y(n-1)
	deallocate(x)
      end