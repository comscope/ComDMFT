      real*8 function f100(fac,j1,j2,j3,m1,m2,m3)
      implicit none
	integer, intent(in) :: j1,j2,j3,m1,m2,m3
      real*8,intent(in) :: fac(*)
	integer :: k1,k2,k3,k4,k5,k6,n1,n2,n,m
	real*8 :: t,t1
      if(m3.ne.m1+m2) goto 2
      k1=j1+j2-j3+1
      k2=j3+j1-j2+1
      k3=j3+j2-j1+1
      k4=j1+j2+j3+2
      t=dfloat(2*j3+1)*fac(k1)*fac(k2)*fac(k3)/fac(k4)
      k1=j1+m1+1
      k2=j1-m1+1
      k3=j2+m2+1
      k4=j2-m2+1
      k5=j3+m3+1
      k6=j3-m3+1
      t=dsqrt(t*fac(k1)*fac(k2)*fac(k3)*fac(k4)*fac(k5)*fac(k6))
      n1=max0(j2-j3-m1,j1-j3+m2,0)+1
      n2=min0(j1+j2-j3,j1-m1,j2+m2)+1
      if(n1.gt.n2) goto 2
      t1=0.d0
      do 1 m=n1,n2
      n=m-1
      k1=j1+j2-j3-n+1
      k2=j1-m1-n+1
      k3=j2+m2-n+1
      k4=j3-j2+m1+n+1
      k5=j3-j1-m2+n+1
  1   t1=t1+dfloat(1+4*(n/2)-2*n)/(fac(m)*fac(k1)*fac(k2)*fac(k3)
     .  *fac(k4)*fac(k5))
      f100=t*t1
      return
  2   f100=0.d0
      return
      end