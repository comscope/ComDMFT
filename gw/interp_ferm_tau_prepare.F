      subroutine interp_ferm_tau_prepare(f_tau,coef,n)
c ------- Prepare for Chebyshev's interpolation ----------
	use manager_mod
	use parallel_mod
	use solid_mod
      implicit none
      integer, intent(in) :: n
	real*8, intent(in) :: f_tau(n,0:n_tau)
	real*8, intent(out) :: coef(n,0:n_tau/2-1,2)
	integer :: i_tau
	real*8, allocatable :: a(:,:)
	allocate(a(n,n_tau/2))
c ---------- Odd part -----------------------------------------
      do i_tau=1,n_tau/2-1
        a(:,i_tau)=f_tau(:,i_tau)-f_tau(:,n_tau-i_tau)
      enddo	
	call dgemm('n','t',n,n_tau/2-1,n_tau/2-1,1.d0,a,n,
     &	         cheb_tau_odd(0,1),n_tau/2-1,0.d0,coef(1,0,1),n)
c ---------- Even part -----------------------------------------
      do i_tau=1,n_tau/2
        a(:,i_tau)=f_tau(:,i_tau)+f_tau(:,n_tau-i_tau)
      enddo	
	call dgemm('n','t',n,n_tau/2,n_tau/2,1.d0,a,n,
     &	         cheb_tau_even(0,1),n_tau/2,0.d0,coef(1,0,2),n)
	deallocate(a)
      end