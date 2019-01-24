      subroutine eig_val_solver_r(n,n0,a,da)
      implicit none
	integer, intent(in) :: n,n0
	real*8, intent(inout) :: a(n0,n0)
	real*8, intent(out) :: da(n)
      integer :: i,lwork,info,j,new
      real*8 :: vl,vu,abstol
      integer, allocatable :: iwork(:),ifail(:)
      real*8, allocatable :: z(:,:),wrk(:)
      allocate(z(n,n)) 
	abstol=1.d-70
      vl=-1.d10
      vu=1.d10
      lwork=max(1,8*n)
      allocate(wrk(lwork),iwork(5*n),ifail(n))
	call DSYEVX('V','A','U',n,a,n0,vl,vu,1,1,abstol,new,da,z,n,wrk,
     &	          lwork,iwork,ifail,info) 
      do j=1,n
        do i=1,n
          a(i,j)=z(i,j)
        enddo
      enddo
	deallocate(wrk,iwork,ifail,z)
      end