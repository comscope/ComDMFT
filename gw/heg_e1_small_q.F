      double precision function heg_e1_small_q(nu,q,ro)
	use units_mod
      implicit none
	real*8, intent(in) :: nu,q,ro
	real*8 :: a
      a=nu*nu+q**4
      heg_e1_small_q=a/(a+16.d0*pi*ro)
      end