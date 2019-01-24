      real*8 function sum_ff(e1,e2)
c     Calculates 1/b Sum_w f(w,e1)f(w,e2)
	use manager_mod
	use solid_mod
      implicit none
	real*8, intent(in) :: e1,e2
	real*8 :: green0
	if(abs(e1-e2)>1.d-8) then
	  sum_ff=(green0(e1,0.d0)-green0(e2,0.d0))/(e1-e2)
	else
	  sum_ff=-betta_t*green0(e1,0.d0)*green0(e1,betta_t)
	endif
      end