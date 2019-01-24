      subroutine dlmtor(r,a,lmax,f,fbar)
c  radial part of damped lmtos f and fbar, l=0 to lmax
	use units_mod
      implicit none
	integer :: l,lmax
      real*8 :: f(0:lmax),fbar(0:lmax),r,a,fbsrpi,z,emz2,erfc0,erfc1,
     &	      erfc2,ta2r,g,flm2,erfc_loc
      fbsrpi=1.d0/dsqrt(pi)
      z=a*r
      emz2=dexp(-z*z)
      erfc0=erfc_loc(z)
      erfc1=-z*erfc0+fbsrpi*emz2
      erfc2=-0.5d0*z*erfc1+0.25d0*erfc0
      f(0)=erfc0/r
      fbar(0)=-erfc2/(a*a*r)
      ta2r=2.d0*a*a*r
      g=2.d0*a*emz2*fbsrpi/r
      flm2=fbsrpi*emz2/z-erfc0
      do l=1,lmax
        f(l)=((l+l-1)/r)*f(l-1)+g
        fbar(l)=((l+l-1)/r)*fbar(l-1)-flm2
        flm2=f(l-1)
        g=g*ta2r
      enddo
      end