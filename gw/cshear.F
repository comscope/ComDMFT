      subroutine cshear(b0,b,ibtr)
c- tries to make microcell more compart by shearing.
c ----------------------------------------------------------------
ci inputs
ci
co outputs
co
cr remarks
cr   ibtr gives the transformation from bo to b.
c ----------------------------------------------------------------
c      implicit none
c passed parameters
c local parameters
      implicit none
	integer :: ibtr(3,3),i,j
      real*8 :: b(3,3),b0(3,3)
      do  8  i = 1, 3
        do  7  j = 1, 3
          b(j,i) = b0(j,i)
          ibtr(j,i) = 0
    7   continue
        ibtr(i,i) = 1
    8 continue
      end
      
      subroutine cross(a,b,c)
**********************************************************
*    cross product (ax,ay,az)=(bx,by,bz)*(cx,cy,cz)      *
**********************************************************
      implicit none
      real*8 a(3),b(3),c(3)
      a(1)=b(2)*c(3)-b(3)*c(2)
      a(2)=b(3)*c(1)-b(1)*c(3)
      a(3)=b(1)*c(2)-b(2)*c(1)
      end
      
      
      double precision function cubint(x0,x1,f0,f1,d0,d1,eps)
      implicit none
      real*8, intent(in) :: x0,x1,f0,f1,d0,d1,eps
	real*8 :: h,d,c,xx,yy,a,b,diskr,b2,a6
      if(x1-x0.lt.eps) then
	  cubint=0.5d0*(x0+x1)
	  return
	endif
      h=x1-x0
      d=f0
      c=d0
      xx=h*(d1-c)
      yy=f1-c*h-d
      a=(xx-2.d0*yy)/h**3
      b=(3.d0*yy-xx)/h/h
	diskr=4.d0*(b*b-3.d0*a*c)
	if(diskr<1.d-30) diskr=1.d-30
      diskr=sqrt(diskr)
      b2=-b-b
      a6=6.d0*a
      cubint=x0+(b2-diskr)/a6
      if(cubint.ge.x0.and.cubint.le.x1) return
      cubint=x0+(b2+diskr)/a6
      end
      
      real*8 function det3(a,b,c)
      implicit none
      real*8 :: a(3),b(3),c(3)
      det3=a(1)*b(2)*c(3)+a(2)*b(3)*c(1)+a(3)*b(1)*c(2)-a(3)*b(2)*c(1)-
     &     a(2)*b(1)*c(3)-a(1)*b(3)*c(2)
      end