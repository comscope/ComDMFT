      real*8 function sum_b0b0(kp1,kp2)
c     Calculates 1/b Sum_v b0(v,kp1)*b0(v,kp2)
	use manager_mod
	use solid_mod
      implicit none
	real*8, intent(in) :: kp1,kp2
	real*8 :: bose_red,bose_red_0
	if(abs(kp1-kp2)>1.d-8) then
	  sum_b0b0=(bose_red(kp1)-bose_red(kp2))/(kp2-kp1)
	else
	  sum_b0b0=bose_red_0(kp1)
	endif
      end