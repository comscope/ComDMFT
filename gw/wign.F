      real*8 function wign(l,m1,m2,al,bet,gam,fac)
      implicit none
      integer :: l,m1,m2
      real*8 :: fac(0:*),al,bet,gam,a,b,d1,d2,a1,a2,dm
      wign=0.0d0
      if(m1.eq.0.and.m2.eq.0) then
         wign=dm(l,0,0,bet,fac)
      else if(m1.eq.0) then
         a=dsqrt(2.d0)*dm(l,0,m2,bet,fac)
         b=m2*gam
         if(m2.lt.0) wign=a*dsin(b)
         if(m2.gt.0) wign=(-1)**m2*a*dcos(b)
      else if(m2.eq.0) then
         a=dsqrt(2.d0)*dm(l,m1,0,bet,fac)
         b=m1*al
         if(m1.lt.0) wign=-a*dsin(b)
         if(m1.gt.0) wign=(-1)**m1*a*dcos(b)
      else
         d1=dm(l,m1,m2,bet,fac)
         d2=dm(l,m1,-m2,bet,fac)
         a=m1*al
         b=m2*gam
         a1=a+b
         a2=a-b
         a=dcos(a1)
         b=dcos(a2)
         a1=dsin(a1)
         a2=dsin(a2)
         if(m1.lt.0) then
            if(m2.lt.0) then
               wign=d1*a-(-1)**m2*d2*b
            else
               wign=-(-1)**m2*(d1*a1+(-1)**m2*d2*a2)
            end if
         else
            if(m2.lt.0) then
               wign=(-1)**m1*(d1*a1-(-1)**m2*d2*a2)
            else
               wign=(-1)**m1*(d1*a*(-1)**m2+d2*b)
            end if
         end if
      end if
      end
