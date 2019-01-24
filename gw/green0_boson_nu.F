      double precision function green0_boson_nu(kap,nu)
c     B(v)=1/(v^2+kap^2)    
	use manager_mod
      implicit none
	real*8, intent(in) :: kap,nu
	green0_boson_nu=1.d0/(nu*nu+kap*kap)
      end