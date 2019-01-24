      complex*16 function sum_fb0(e,kp,om)
c     Calculates 1/b Sum_v f(w+v,e)b0(v,kp)
	use manager_mod
	use solid_mod
      implicit none
	real*8, intent(in) :: e,kp,om
	real*8 :: bose_red,green0
	complex*16 :: fermi_omega
	sum_fb0=fermi_omega(e+kp,om)*(bose_red(kp)-green0(e,0.d0)-0.5d0)
      end