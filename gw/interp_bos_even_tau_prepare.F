      subroutine interp_bos_even_tau_prepare(f_tau,coef,n)
c ------- Prepare for Chebyshev's interpolation ----------
	use manager_mod
	use parallel_mod
	use solid_mod
      implicit none
      integer, intent(in) :: n
	real*8, intent(in) :: f_tau(n,0:n_tau/2)
	real*8, intent(out) :: coef(n,0:n_tau/2-1)
	call dgemm('n','t',n,n_tau/2,n_tau/2,1.d0,f_tau(1,1),n,
     &	         cheb_tau_even(0,1),n_tau/2,0.d0,coef(1,0),n)
      end