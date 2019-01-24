      SUBROUTINE SGAUNT
************************************************************
*   Calculation of the Gaunt coefficients C(L2M2,L1M1,LM)  *
*                                                          *
*    l2m2                    /         *                   *
*   C       =C(l2m2,l1m1,lm)=\dr*Y(r)*Y(r)*Y(r)            *
*    lm,l1m1                 /    lm   l1m1 l2m2           *
*                                                          *
*    and C.ne.0 when l2=/l1-l/,/l1-l/+2,...,l1+l,m2=m1-m   *
*    Y(lm) is a complex spherical garmonic with a phase    *
*    after Condon and Shortley                             *
* Written by S.Yu.Savrasov (P.N.Lebedev Physical Institute)*
************************************************************
	use atom_mod
	use manager_mod
	use units_mod
      IMPLICIT none
	integer :: l1,m1,l1m1,l,m,lm,l2,ll2,m2
	real*8 :: aj,bj,am,bm,cj,cm,clebsch,dl1,dl2,dl3,a1,a2
      allocate(gaunt(0:maxpb/2,(maxb+2)**2,(maxb+2)**2))
      DO L1=0,maxb+1
        DO M1=-L1,L1
          L1M1=L1*(L1+1)+M1+1
          DO L =0,maxb+1
            DO M =-L,L
              LM=L*(L+1)+M+1
              DO L2=ABS(L1-L),min(L1+L,maxpb),2
                LL2=L2/2
                M2=M1-M
                IF(ABS(M2).LE.L2)THEN     !!! selection rule
                  AJ=L
                  BJ=L2
                  AM=M
                  BM=M2
                  CJ=L1
                  CM=M1
                  A1=CLEBSCH(AJ,BJ,0.D0,0.D0,CJ,0.D0)  
                  A2=CLEBSCH(AJ,BJ,AM,BM,CJ,CM)        
                  DL1=2*L +1
                  DL2=2*L2+1
                  DL3=2*L1+1
                  gaunt(LL2,L1M1,LM)=A1*A2*SQRT(DL1*DL2/DL3/4.D0/PI)
                ELSE IF(ABS(M2).GT.L2)THEN
                  gaunt(LL2,L1M1,LM)=0.d0
                ENDIF
	        enddo  !! over l2
	      enddo  !! over m
	    enddo  !! over l
	  enddo !! over m1
	enddo   !! l1
      END