      subroutine spline_hmg(h,y,a,b,c,d,n,key,dy_0,dy_n)
c	key = 0 - Free boundary conditions (y''(0)=y''(n)=0)
c	key = 1 - y'(0) and y'(n) are given
c	key = 2 - y''(0)=0 and y'(n) is given
c	key = 3 - y'(0) is given and y''(n)=0
      implicit none
	integer, intent(in) :: n,key
	real*8, intent(in) :: h,y(0:n),dy_0,dy_n
	real*8, intent(out) :: a(n),b(n),c(n),d(n)
	real*8, allocatable :: x(:)
	integer :: i
c -------------------------------------------------
	a(1)=0.d0
	a(2:n-1)=h/3.d0
c -------------------------------------------------
	do i=1,n-1
	  c(i)=h/3.d0
	enddo
	c(n-1)=0.d0
c -------------------------------------------------
	do i=2,n-2
	  b(i)=-4.d0*h/3.d0
	  d(i)=(y(i+1)-y(i))/h-(y(i)-y(i-1))/h
	enddo
c -------------------------------------------------
	if(key==0.or.key==2) then
	  b(1)=-4.d0*h/3.d0
	  d(1)=(y(2)-y(1))/h-(y(1)-y(0))/h
	else if(key==1.or.key==3) then
	  b(1)=-h/2.d0-2.d0*h/3.d0
	  d(1)=-1.5d0*(y(1)-y(0))/h+0.5d0*dy_0+(y(2)-y(1))/h
	endif
c -------------------------------------------------
	if(key==0.or.key==3) then
	  b(n-1)=-4.d0*h/3.d0
	  d(n-1)=(y(n)-y(n-1))/h-(y(n-1)-y(n-2))/h
	else if(key==1.or.key==2) then
	  b(n-1)=-h/2.d0-2.d0*h/3.d0
	  d(n-1)=1.5d0*(y(n)-y(n-1))/h-0.5d0*dy_n
     &	                             -(y(n-1)-y(n-2))/h
	endif
	allocate(x(n))
	call spline(a,b,c,d,n-1,x(2))
	c(2:)=x(2:)
	if(key==0.or.key==2) then
	  c(1)=0.d0
	else if(key==1.or.key==3) then
	  c(1)=1.5d0*(y(1)-y(0))/h/h-1.5d0*dy_0/h-0.5d0*c(2)
	endif
	do i=1,n-1
	  b(i)=(y(i)-y(i-1))/h-h*(c(i+1)+2.d0*c(i))/3.d0
	  d(i)=(c(i+1)-c(i))/3.d0/h
	  a(i)=y(i-1)
	enddo
	a(n)=y(n-1)
	if(key==0.or.key==3) then
	  b(n)=(y(n)-y(n-1))/h-2.d0*h*c(n)/3.d0
	  d(n)=-c(n)/3.d0/h
	else if(key==1.or.key==2) then
	  b(n)=1.5d0*(y(n)-y(n-1))/h-0.5d0*dy_n-0.5d0*h*c(n)
	  d(n)=(dy_n-b(n)-2.d0*h*c(n))/3.d0/h/h
	endif
	deallocate(x)
      end