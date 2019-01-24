      complex*16 function fermi_omega(e,w)
c     F(w)=1/(iw-e)
	use manager_mod
	use solid_mod
      implicit none
	real*8, intent(in) :: w,e
	fermi_omega=(1.d0,0.d0)/dcmplx(-e,w)
      end