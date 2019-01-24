      subroutine x_integral_sigma_0(a,b,q,k,wv,chem)
c ------- Calculates X-integral for noninteracting G -----------
      use heg_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	real*8, intent(in) :: k,q,wv,chem
	real*8, intent(out) :: a,b
	real*8 :: kq,d,c1,c2,den
	d=wv*wv
	if(k<1.d-10) then
	  c1=chem-q*q
	  den=d+c1*c1
	  a=4.d0*c1/den
	  b=-4.d0*wv/den
	else
	  kq=k*q
	  c1=chem-(k+q)**2
	  c2=chem-(k-q)**2
	  den=d+c1*c2
	  if(8.d0*kq*abs(c1+2.d0*kq)/d<1.d-7) then
	    a=4.d0*(c1+2.d0*kq)/(d+c2*c2)
	  else
	    a=-0.5d0/kq*log((d+c1*c1)/(d+c2*c2))
	  endif
	  b=atan(-4.d0*kq*abs(wv)/den)
	  if(chem>0.d0) b=b+pi/2.d0*(1.d0-sign(1.d0,den))*sign(1.d0,c1)
	  b=sign(1.d0,wv)*b/kq
	endif
      end