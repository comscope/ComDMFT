      subroutine spline(a,b,c,d,n,x)
      implicit none
	integer, intent(in) :: n
	real*8, intent(in) :: a(n),b(n),c(n),d(n)
	real*8, intent(out) :: x(n)
	real*8, allocatable :: ksi(:),eta(:)
	integer :: i
	real*8 :: tmp
	allocate(ksi(n+1),eta(n+1))
c -------- Direct running -----------------------------------------------
	ksi(1)=0.d0
	eta(1)=0.d0
	do i=1,n
	  tmp=b(i)-a(i)*ksi(i)
	  ksi(i+1)=c(i)/tmp
	  eta(i+1)=(a(i)*eta(i)-d(i))/tmp
	enddo
c -------- Backward running -----------------------------------------------
      x(n)=eta(n+1)
	do i=n-1,1,-1
	  x(i)=ksi(i+1)*x(i+1)+eta(i+1)
	enddo
	deallocate(ksi,eta)
      end
      subroutine spline_c(a,b,c,d,n,x)
      implicit none
	integer, intent(in) :: n
	complex*16, intent(in) :: a(n),b(n),c(n),d(n)
	complex*16, intent(out) :: x(n)
	complex*16, allocatable :: ksi(:),eta(:)
	integer :: i
	complex*16 :: tmp
	allocate(ksi(n+1),eta(n+1))
c -------- Direct running -----------------------------------------------
	ksi(1)=(0.d0,0.d0)
	eta(1)=(0.d0,0.d0)
	do i=1,n
	  tmp=b(i)-a(i)*ksi(i)
	  ksi(i+1)=c(i)/tmp
	  eta(i+1)=(a(i)*eta(i)-d(i))/tmp
	enddo
c -------- Backward running -----------------------------------------------
      x(n)=eta(n+1)
	do i=n-1,1,-1
	  x(i)=ksi(i+1)*x(i+1)+eta(i+1)
	enddo
	deallocate(ksi,eta)
      end