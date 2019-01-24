      real*8 function integral_pw_i(q,l,r,bess)
      integer, intent(in) :: l
      real*8, intent(in) :: q,r,bess(0:l+1)
	if(q<1.d-10) then
	  if(l==0) integral_pw_i=r**3/3.d0
	  if(l/=0) integral_pw_i=0.d0
	else
	  integral_pw_i=r**(l+2)*bess(l+1)/q
      endif
      end