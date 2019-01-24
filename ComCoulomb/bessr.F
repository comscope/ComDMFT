      SUBROUTINE BESSR(k,r,BES,BESD,LMAX)
******************************************************************
*     Caclculate Bessel functions of a real argument               *
*     as well as their derivatives with respect to argument.       *
******************************************************************
      IMPLICIT none
      integer :: lmax,li,il,l,l00,lup
      double precision ::
     $  k,x,BESW(0:26),BES(0:LMAX),BESD(0:LMAX),ww,bnor,sqbn,r,
     &  fac2
C     
c     IF(LMAX.GT.26)STOP'<BESSC> LMAX>26 DOES NOT WORK'
c     IF(LMAX.GT.20)PRINT*,'LMAX EXCEEDED 20 IN BESSR FUNCTION'
C     
      x=k*r
      IF (ABS(X).GT.10.D0) THEN ! DIRECT RECURSION
        BES(0)=SIN(X)/X
        BESD(0)=-(SIN(X)/X-COS(X))/X
        IF(LMAX.GT.0)BES(1)=-BESD(0)
        IF(LMAX.GT.0)BESD(1)=BES(0)-2.D0*BES(1)/X
        DO LI=0,LMAX-2
          BES(LI+2)=DBLE(2*LI+3)/X*BES(LI+1)-BES(LI)
          BESD(LI+2)=BES(LI+1)-DBLE(LI+3)*BES(LI+2)/X
        enddo
      ELSE                      ! REVERSE RECURSION
        DO IL=0,LMAX
          BES(IL)=0.D0
          BESD(IL)=0.D0
        enddo
        IF(ABS(X).LT.1.D-08) THEN
          BES(0)=1.D0
          BESD(0)=0.D0
          do l=1,lmax
            if(l.eq.1) then
              besd(l)=k/3.d0
            else
              besd(l)=l*k**l*r**(l-1)/fac2(l+l+1)
            endif
            bes(l)=r*besd(l)/l
          enddo
          RETURN
        ENDIF
        IF(ABS(X).EQ.1.D0) X=1.D0+1.D-06
        L00=-24.D0/LOG10(ABS(X))
        LUP=MIN(26,L00)
        IF (L00.LE.0) LUP=26
 1      WW=X**LUP/FAC2(2*LUP+1)*(1.D0-X**2/2.D0/DBLE(2*LUP+3))
        IF(ABS(WW).EQ.0.D0) THEN
          LUP=LUP-1
          GOTO 1
        ENDIF
        BESW(LUP)=WW
        BESW(LUP-1)=X**(LUP-1)/FAC2(2*LUP-1)*(1.D0-
     &    X**2/2.D0/DBLE(2*LUP+1))
        BNOR=DBLE(2*LUP+1)*BESW(LUP)*BESW(LUP)+
     &    DBLE(2*LUP-1)*BESW(LUP-1)*BESW(LUP-1)
        DO LI=LUP-2,0,-1
          BESW(LI)=DBLE(2*LI+3)/X*BESW(LI+1)-BESW(LI+2)
          BNOR=BNOR+DBLE(2*LI+1)*BESW(LI)*BESW(LI)
        enddo
        SQBN=SQRT(BNOR)
        BES(0)=BESW(0)/SQBN
        BESD(0)=(COS(X)-SIN(X)/X)/X
        DO LI=1,MIN(LMAX,LUP)
          BES(LI)=BESW(LI)/SQBN
          BESD(LI)=BES(LI-1)-DBLE(LI+1)*BES(LI)/X
        enddo
      ENDIF
      do l=0,lmax
        besd(l)=k*besd(l)
      enddo
      END
