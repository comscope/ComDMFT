	SUBROUTINE intalt(N,f,g)
      IMPLICIT none
	integer :: n,i
      real*8 :: F(n),g(n),h12,h3,sum,f1,f2,f3
C     N - NUMBER OF points ( odd !!! )
C     g(I) = INTEGRAL FROM 0 TO I F(Y)dY
      H12=1.d0/12.D0
      H3=1.d0/3.D0
      SUM=0.D0
      g(1)=0.D0
      DO I=1,N-2,2
        F1=F(i)
        F2=F(I+1)
        F3=F(I+2)
        g(I+1)=SUM+H12*(5.D0*F1+8.D0*F2-F3)
        g(I+2)=SUM+H3*(F1+4.D0*F2+F3)
        SUM=g(I+2)
      ENDDO
      END