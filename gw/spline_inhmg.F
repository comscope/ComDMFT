      subroutine spline_inhmg(mesh,y,a,b,c,d,n,key,dy_0,dy_n)
c	key = 0 - Free boundary conditions (y''(0)=y''(n)=0)
c	key = 1 - y'(0) and y'(n) are given
c	key = 2 - y''(0)=0 and y'(n) is given
c	key = 3 - y'(0) is given and y''(n)=0
      implicit none
	integer, intent(in) :: n,key
	real*8, intent(in) :: mesh(0:n),y(0:n),dy_0,dy_n
	real*8, intent(out) :: a(n),b(n),c(n),d(n)
	real*8, allocatable :: x(:),h(:)
	integer :: i
	allocate(h(n))
	do i=1,n
	  h(i)=mesh(i)-mesh(i-1)
	enddo
c -------------------------------------------------
	a(1)=0.d0
	a(2:n-1)=h(2:n-1)/3.d0
c -------------------------------------------------
	do i=1,n-1
	  c(i)=h(i+1)/3.d0
	enddo
	c(n-1)=0.d0
c -------------------------------------------------
	do i=2,n-2
	  b(i)=-2.d0*(h(i)+h(i+1))/3.d0
	  d(i)=(y(i+1)-y(i))/h(i+1)-(y(i)-y(i-1))/h(i)
	enddo
c -------------------------------------------------
	if(key==0.or.key==2) then
	  b(1)=-2.d0*(h(1)+h(2))/3.d0
	  d(1)=(y(2)-y(1))/h(2)-(y(1)-y(0))/h(1)
	else if(key==1.or.key==3) then
	  b(1)=-h(1)/2.d0-2.d0*h(2)/3.d0
	  d(1)=-1.5d0*(y(1)-y(0))/h(1)+0.5d0*dy_0+(y(2)-y(1))/h(2)
	endif
c -------------------------------------------------
	if(key==0.or.key==3) then
	  b(n-1)=-2.d0*(h(n-1)+h(n))/3.d0
	  d(n-1)=(y(n)-y(n-1))/h(n)-(y(n-1)-y(n-2))/h(n-1)
	else if(key==1.or.key==2) then
	  b(n-1)=-h(n)/2.d0-2.d0*h(n-1)/3.d0
	  d(n-1)=1.5d0*(y(n)-y(n-1))/h(n)-0.5d0*dy_n
     &	                             -(y(n-1)-y(n-2))/h(n-1)
	endif
	allocate(x(n))
	call spline(a,b,c,d,n-1,x(2))
	c(2:)=x(2:)
	if(key==0.or.key==2) then
	  c(1)=0.d0
	else if(key==1.or.key==3) then
	  c(1)=1.5d0*(y(1)-y(0))/h(1)/h(1)-1.5d0*dy_0/h(1)-0.5d0*c(2)
	endif
	do i=1,n-1
	  b(i)=(y(i)-y(i-1))/h(i)-h(i)*(c(i+1)+2.d0*c(i))/3.d0
	  d(i)=(c(i+1)-c(i))/3.d0/h(i)
	  a(i)=y(i-1)
	enddo
	a(n)=y(n-1)
	if(key==0.or.key==3) then
	  b(n)=(y(n)-y(n-1))/h(n)-2.d0*h(n)*c(n)/3.d0
	  d(n)=-c(n)/3.d0/h(n)
	else if(key==1.or.key==2) then
	  b(n)=1.5d0*(y(n)-y(n-1))/h(n)-0.5d0*dy_n-0.5d0*h(n)*c(n)
	  d(n)=(dy_n-b(n)-2.d0*h(n)*c(n))/3.d0/h(n)/h(n)
	endif
	deallocate(x,h)
      end