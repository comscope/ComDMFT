      real*8 function green_34_tau(a,t)
c     Calculates correlation part of Green's function of the form
c       {(b-t)^2/2-b(t+5b/2)G(a,b)+b^2G^2(a,b)}G(a,b-t)
	use manager_mod
      implicit none
	real*8, intent(in) :: a,t
	real*8 :: green0,gb,gbt,b
	b=betta_t
	gb=green0(a,b)
	gbt=green0(a,b-t)
	green_34_tau=((b-t)**2/2.d0-b*(t+2.5d0*b)*gb+(b*gb)**2)*gbt
      end