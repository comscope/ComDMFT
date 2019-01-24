      real*8 function integral_pw_j(q,l,r,sm,bess,bess1)
      integer, intent(in) :: l
      real*8, intent(in) :: q,r,sm,bess(0:l),bess1(0:l)
	if(q<1.d-10) then
	  if(l==0) integral_pw_j=(sm*sm-r*r)/2.d0
	  if(l/=0) integral_pw_j=0.d0
	else
	  if(l==0) integral_pw_j=(cos(q*r)-cos(q*sm))/q**2
	  if(l/=0) integral_pw_j=(bess(l-1)/r**(l-1)
     &	                       -bess1(l-1)/sm**(l-1))/q
      endif
      end