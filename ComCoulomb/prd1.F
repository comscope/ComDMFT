      subroutine prd1(a,b,c)
      integer :: i2,i3
      double precision :: a(3),b(3,3),c(3)
      do 20 i2=1,3
        c(i2)=0.d0
        do 10 i3=1,3
          c(i2)=c(i2)+a(i3)*b(i3,i2)
 10     continue
 20   continue
      end
