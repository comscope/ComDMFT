      real*8 FUNCTION CLEBSCH(AJ,BJ,AM,BM,CJ,CM)
******************************************************************
* Program calculates Clebsch-Gordan coefficients                 *
* See: Landau and Lifshitz, Vol.3                                *
*      cj,cm                                                     *
*     C                                                          *
*      aj,am,bj,bm                                               *
* Written by A.Soldatov (IAE)                                    *
******************************************************************
      IMPLICIT none
	real*8, intent(in) :: aj,bj,am,bm,cj,cm
	integer :: n,k,i,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10
      real*8 :: F(100),x,s,e,c
      SAVE F,N,K,X
      DATA N/100/, K/0/,  X/2.D0/
      F(1)=0.D0
      F(2)=0.D0
	if(k<=0) then
	  k=1
        DO I=3,N
          F(I)=F(I-1)+DLOG(X)
          X=X+1.D0
	  enddo
	endif
      I=AM+BM-CM+.1D0
      IF(I/=0) goto 100
      I1=AJ+BJ-CJ+1.1D0
      IF(I1<=0) goto 100
      I2=AJ-BJ+CJ+1.1D0
      IF(I2<=0) goto 100
      I3=BJ+CJ-AJ+1.1D0
      IF(I3<=0) goto 100
      X=AJ+BJ+CJ+2.1D0
      I4=X
      I=X+.6D0
      I=I4-I
      IF(I/=0) goto 100
      X=AJ+AM+1.1D0
      I5=X
      IF(I5<=0) goto 100
      I=X+.6D0
      I=I-I5
      IF(I/=0) goto 100
      I6=AJ-AM+1.1D0
      IF(I6<=0) goto 100
      X=BJ+BM+1.1D0
      I7=X
      IF(I7<=0) goto 100
      I=X+.6D0
      I=I-I7
      IF(I/=0) goto 100
      I8=BJ-BM+1.1D0
      IF(I8<=0) goto 100
      X=CJ+CM+1.1D0
      I9=X
      IF(I9<=0) goto 100
      I=X+.6D0
      I=I-I9
      IF(I/=0) goto 100
      I10=CJ-CM+1.1D0
      IF(I10<=0) goto 100
      X=F(I1)+F(I2)+F(I3)-F(I4)
      I=I5-I6
      IF(I/=0) goto 18
      I=I7-I8
      IF(I==0) goto 200
   18 X=X+F(I5)+F(I6)+F(I7)+F(I8)+F(I9)+F(I10)
      X=X*.5D0
      I10=MIN0(I1,I6,I7)
      I2=I1-I5
      I3=I1-I8
      I9=MAX0(0,I2,I3)+1
      I1=I1+1
      I6=I6+1
      I7=I7+1
      I8=I9/2
      E=1.D0
      I5=I9*.5D0+.6D0
      I8=I8-I5
      IF(I8==0) E=-1.D0
      S=0.D0
      DO I=I9,I10
        C=X-F(I)-F(I1-I)-F(I6-I)-F(I7-I)-F(I-I2)-F(I-I3)
        S=S+E*DEXP(C)
        E=1.D0-E-1.D0
      enddo
      CLEBSCH=DSQRT(CJ+CJ+1.D0)*S
      RETURN
  200 I=I4/2
      I5=I4*.5D0+.6D0
      I=I-I5
      IF(I/=0) goto 100
      I6=I5-I6+1
      I7=I5-I8+1
      I8=I5-I10+1
      S=X*0.5D0+F(I5)-F(I6)-F(I7)-F(I8)
      S=DEXP(S)
      I5=I8/2
      I6=I8*.5D0+.6D0
      I5=I5-I6
      IF(I5==0) S=1.D0-S-1.D0
      CLEBSCH=S*DSQRT(CJ+CJ+1.D0)
      RETURN
  100 CLEBSCH=0.D0
      RETURN
      END