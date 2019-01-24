      subroutine boson_cheb_coeff_nu(a_nu,a_cheb,n)
c ------- Chebyshev's coefficients for A (bosonic) --------------------
c --- Input and output arrays are in compact form (triangle if needed)
c --- where n - total number of elements ------------------------------
	use manager_mod
	use parallel_mod
	use solid_mod
      implicit none
	integer, intent(in) :: n
	real*8, intent(in) :: a_nu(n,0:n_nu)
	real*8, intent(out) :: a_cheb(n,0:n_nu_asy-2)
	integer :: ic_nu,l,i_nu
	a_cheb=0.d0
	do ic_nu=1,n_nu_asy-1
	  i_nu=n_nu_exa+ic_nu
	  do l=0,n_nu_asy-2
	    a_cheb(:,l)=a_cheb(:,l)+cheb_nu(l,ic_nu)*a_nu(:,i_nu)
	  enddo
	enddo
      end