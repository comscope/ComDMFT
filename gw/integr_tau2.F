      subroutine integr_tau2(ff,dd,i_tau,n,tran)
c ------ Integral of f(t';t) where t' in t ---------
	use manager_mod
	use parallel_mod
	use solid_mod
      implicit none
      integer, intent(in) :: i_tau,n
      real*8, intent(in) :: ff(n,0:n_tau),tran(n_tau/2)
      real*8, intent(out) :: dd(n)
	integer :: m,m1,mm,nn,i
	real*8, allocatable :: f(:,:)
	m=n_bound_tau(i_tau)
      m1=n_tau1(i_tau)-n_bound_tau(i_tau)
      mm=n_tau1(i_tau)
      nn=m/2+m1/2
	allocate(f(n,nn))
	do i=1,m/2
	  f(:,i)=ff(:,i)+ff(:,m-i)
	enddo
	do i=1,m1/2
	  f(:,m/2+i)=ff(:,m+i)+ff(:,mm-i)
	enddo
	call dgemm('n','n',n,1,nn,1.d0,f,n,tran,n_tau/2,0.d0,dd,n)
      deallocate(f)
      end