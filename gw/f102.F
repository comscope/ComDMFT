      real*8 function f102(fac,l1,l2,l3)
      implicit none
	integer, intent(in) :: l1,l2,l3
      real*8, intent(in) :: fac(*)
      integer :: x,p,lt
      lt=l1+l2+l3
      p=lt/2
      if(2*p.ne.lt) goto 1
      f102=dsqrt(dfloat(2*l3+1)/dfloat(lt+1))
      f102=f102*fac(p+1)/dsqrt(fac(2*p+1))
      x=p-l1
      f102=f102*dsqrt(fac(2*x+1))/fac(x+1)
      x=p-l2
      f102=f102*dsqrt(fac(2*x+1))/fac(x+1)
      x=p-l3
      f102=f102*dsqrt(fac(2*x+1))/fac(x+1)
      if(x.gt.2*(x/2)) f102=-f102
      return
  1   f102=0.d0
      return
      end