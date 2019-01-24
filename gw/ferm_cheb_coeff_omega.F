      subroutine ferm_cheb_coeff_omega(ab_omega,ab_cheb,n)
c ------- Chebyshev's coefficients for AB (fermionic) -----------------
c --- Input and output arrays are in compact R*8 form -----------------
	use manager_mod
	use parallel_mod
	use solid_mod
      implicit none
	integer, intent(in) :: n
	real*8, intent(in) :: ab_omega(n,n,2,0:n_omega)
	real*8, intent(out) :: ab_cheb(n,n,2,0:n_omega_asy-2)
	integer :: ic_omega,l,i_omega
	ab_cheb=0.d0
	do ic_omega=1,n_omega_asy-1
	  i_omega=n_omega_exa+ic_omega
	  do l=0,n_omega_asy-2
	    ab_cheb(:,:,:,l)=ab_cheb(:,:,:,l)+cheb_omega(l,ic_omega)
     &	                                   *ab_omega(:,:,:,i_omega)
	  enddo
	enddo
      end