      double precision function det3(a,b,c)
      implicit none
      double precision :: a(3),b(3),c(3)
      det3=a(1)*b(2)*c(3)+a(2)*b(3)*c(1)+a(3)*b(1)*c(2)-a(3)*b(2)*c(1)-
     &  a(2)*b(1)*c(3)-a(1)*b(3)*c(2)
      end
