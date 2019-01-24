      complex*16 function green_34_omega(a,omega)
c     Calculates correlation part of Green's function of the form
c                          1/(iw+a)^3      
	use manager_mod
      implicit none
	real*8, intent(in) :: a,omega
	green_34_omega=(1.d0,0.d0)/dcmplx(a,omega)**3
      end