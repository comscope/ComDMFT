      subroutine prd(a,b,c,l)
	integer :: l,m1,m2,m3
      real*8 :: a(-l:l,-l:l),b(-l:l,-l:l),c(-l:l,-l:l)
      do m1=-l,l
        do m2=-l,l
          c(m1,m2)=0.d0
          do m3=-l,l
            c(m1,m2)=c(m1,m2)+a(m1,m3)*b(m3,m2)
	    enddo
	  enddo
	enddo
      end
	subroutine prd_j(a,b,c,j)
	integer :: j,m1,m2,m3
      complex*16 :: a(j+1,j+1),b(j+1,j+1),c(j+1,j+1)
      do m1=1,j+1
        do m2=1,j+1
          c(m1,m2)=(0.d0,0.d0)
          do m3=1,j+1
            c(m1,m2)=c(m1,m2)+a(m1,m3)*b(m3,m2)
	    enddo
	  enddo
	enddo
      end
