      double precision function bose_red(kap)
c     B(kap)=1/2 * (1+e^(-kap*b))/(1-e^(-kap*b))      
	use manager_mod
      implicit none
	real*8, intent(in) :: kap
	real*8 :: a,b,t
	if(kap>=0) then
	  t=kap*betta_t
	  a=1.d0
	  if(t<100.d0) a=1.d0+exp(-t)
	  b=1.d0
	  if(t<100.d0) b=1.d0-exp(-t)
	else ! if(kap<0) then
	  t=kap*betta_t
	  a=1.d0
	  if(abs(t)<100.d0) a=exp(t)+1.d0
	  b=1.d0
	  if(abs(t)<100.d0) b=exp(t)-1.d0
	endif
	bose_red=0.5d0*a/b
      end
      
      
      double precision function bose_red_0(kap)
c     B(kp)=b * e^(-kap*b)/(1-e^(-kap*b))**2      
	use manager_mod
      implicit none
	real*8, intent(in) :: kap
	real*8 :: a,b,t
	if(kap>=0) then
	  t=kap*betta_t
	  a=0.d0
	  if(t<100.d0) a=exp(-t)
	  b=1.d0
	  if(t<100.d0) b=1.d0-exp(-t)
	else ! if(kap<0) then
	  t=kap*betta_t
	  a=0.d0
	  if(abs(t)<100.d0) a=exp(t)
	  b=1.d0
	  if(abs(t)<100.d0) b=exp(t)-1.d0
	endif
	bose_red_0=betta_t*a/b/b
      end
