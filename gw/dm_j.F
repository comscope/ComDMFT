      real*8 function dm_j(j,m1,m2,b,fac)
      implicit none
      integer :: j,m1,m2,k,k1,k2,k3,k4,ma,na,la,i
      real*8 :: fac(0:*),b,bb,sb,cb,d,pja
      bb=0.5d0*b
      sb=sin(bb)
      cb=cos(bb)
      k1=(j+m2)/2
      k2=(j-m2)/2
      k3=(j+m1)/2
      k4=(j-m1)/2
      k=amin0(k1,k2,k3,k4)
      if(k.eq.k1) then
        ma=(m1-m2)/2
        na=(-m1-m2)/2
        la=ma
      else if(k.eq.k2) then
        ma=(m2-m1)/2
        na=(m1+m2)/2
        la=0
      else if(k.eq.k3) then
        ma=(m2-m1)/2
        na=(-m1-m2)/2
        la=0
      else
        ma=(m1-m2)/2
        na=(m1+m2)/2
        la=ma
      end if
      d=dsqrt(fac(k)*fac(j-k)/(fac(k+ma)*fac(k+na)))
      d=(-1)**la*d
      do i=1,ma
        d=d*sb
      enddo
      do i=1,na
        d=d*cb
      enddo
      dm_j=d*pja(k,ma,na,cos(b),fac)
      end
