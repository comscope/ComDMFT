      subroutine eig_val_solver(n,n0,a,da)
      implicit none
	integer, intent(in) :: n,n0
	complex*16, intent(inout) :: a(n0,n0)
	real*8, intent(out) :: da(n)
      integer :: i,lwork,info,j,new
      real*8 :: vl,vu,abstol
      integer, allocatable :: iwork(:),ifail(:)
      real*8, allocatable :: rwork(:)
	complex*16, allocatable :: z(:,:),work(:)
      allocate(z(n,n)) 
      abstol=1.d-70
      vl=-1.d10
      vu=1.d10
      lwork=2*max(1,2*n-1)
      allocate(work(lwork),rwork(7*n),iwork(5*n),ifail(n))
      call zheevx('V','A','L',n,a,n0,vl,vu,1,1,abstol,new,da,z,n,work,
     &            lwork,rwork,iwork,ifail,info)
      do j=1,n
        do i=1,n
          a(i,j)=z(i,j)
        enddo
      enddo
	deallocate(work,rwork,iwork,ifail,z)
      end
