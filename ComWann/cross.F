      subroutine dcross(a,b,c)
**********************************************************
*     cross product (ax,ay,az)=(bx,by,bz)*(cx,cy,cz)      *
**********************************************************
      implicit none
      double precision, intent(in) :: b(3),c(3)
      double precision, intent(out) :: a(3)      
      a(1)=b(2)*c(3)-b(3)*c(2)
      a(2)=b(3)*c(1)-b(1)*c(3)
      a(3)=b(1)*c(2)-b(2)*c(1)
      end



      double precision function dcrossdot(a, b, c)
      implicit none
      double precision, intent(in) :: a(3),b(3),c(3)
      double precision :: d(3)
      
      call dcross(d, b, c)
      dcrossdot=sum(d*a)
      end
