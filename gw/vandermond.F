      subroutine vandermond(a,n,ksi)
      implicit none
      integer, intent(in) :: n
	real*8, intent(in) :: ksi(n)
	real*8, intent(out) :: a(n,n)
	integer :: i,j,info
	integer, allocatable :: ipiv(:)
	real*8, allocatable :: work(:,:)
c ------------------------------------------------------------
      a(:,1)=1.d0
	do j=2,n
	  do i=1,n
	    a(i,j)=a(i,j-1)*ksi(i)
	  enddo
	enddo
	allocate(ipiv(n),work(n,n))
	call dgetrf(n,n,a,n,ipiv,info)
	call dgetri(n,a,n,ipiv,work,n**2,info)
	deallocate(ipiv,work)
      end