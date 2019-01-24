      real*8 function integral_pw_k(q,q1,l,sm,bess,bess1)
      integer, intent(in) :: l
      real*8, intent(in) :: q,q1,sm,bess(0:l+2),bess1(0:l+2)
	if(q1<1.d-10) then
	  if(q<1.d-10) then
	    if(l==0) integral_pw_k=sm**5/7.5d0
	    if(l/=0) integral_pw_k=0.d0
	  else
	    if(l==0) integral_pw_k=sm**3*(q*sm*bess(1)+bess(2))/3.d0/q/q
	    if(l/=0) integral_pw_k=0.d0
	  endif
	else
	  if(q<1.d-10) then
	    if(l==0) integral_pw_k=sm**3*(q1*sm*bess1(1)+bess1(2))/3.d0/q1
     & 	                                                             /q1
	    if(l/=0) integral_pw_k=0.d0
	  else if(abs(q-q1)<1.d-10) then
	    integral_pw_k=sm**3*((2*l+3)*bess(l+1)**2
     &	                      -(2*l+1)*bess(l)*bess(l+2))/2.d0/q/q
        else
	    integral_pw_k=sm**3*(bess(l+1)*bess1(l+1)/q/q1
     &	                      +(2*l+1)*(bess(l+2)*bess1(l)
     &                                 -bess(l)*bess1(l+2))
     &                        /(2.d0*l+3)/(q*q-q1*q1))
	  endif
      endif
      end