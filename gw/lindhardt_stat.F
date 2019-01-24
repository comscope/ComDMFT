      double precision function lindhardt_stat(x)
      implicit none
	real*8, intent(in) :: x
	real*8 :: a
	if(x==0.d0) then
	  lindhardt_stat=1.d0
	else
	  a=(1.d0-x*x)*log(abs((1.d0+x)/(1.d0-x)))
	  lindhardt_stat=0.5d0*(1.d0+0.5d0*a/x)
	endif
      end
      
      
      
      	
	double precision function lindhard(q,u,kf)
	use units_mod
	implicit none
	real*8, intent(in) :: q,u,kf
	real*8 :: u2,a,qp1,qm1,b,c,x
	if(q<1.d-9.and.u<1.d-9) then
	  lindhard=-kf/2.d0/pi/pi
	else if(q<1.d-9) then
	  lindhard=0.d0
	else if(u<1.d-9.and.abs(q-1.d0)<1.d-12) then
	  a=-1.d0
	  lindhard=kf/4.d0/pi/pi*a
	else
	  u2=u*u
	  qp1=1.d0+q
	  qm1=1.d0-q
	  b=u2+qm1**2
	  c=q*q-u2-1.d0
	  if(q<1.d-5.or.u>100.d0) then
	    x=4.d0*q/b
	    a=-(1.d0+2.d0*q*qm1/b)*(1.d0-x/2+x*x/3-x*x*x/4)
	  else
	    a=log((u2+qp1**2)/b)
	    a=c/4.d0/q*a
	  endif
	  a=a-1.d0
	  if(u>=1.d-9) a=a+u*(atan(qp1/u)+atan(qm1/u))
	  lindhard=kf/4.d0/pi/pi*a
	endif
	end