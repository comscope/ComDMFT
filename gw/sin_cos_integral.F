      subroutine sin_cos_integral(x,ci,si)   
      implicit none
      real*8, intent(in) :: x
      real*8, intent(out) :: ci,si
      integer :: k,m
      real*8 :: p2,el,eps,x2,xr,xa1,xa0,xa,xs,xg1,xg2,xcs,xss,xf,xg
      real*8, allocatable :: bj(:)
      allocate(BJ(101)) 
      P2=1.570796326794897D0 
      EL=.5772156649015329D0 
      EPS=1.D-15 
      X2=X*X 
      IF (X.EQ.0.0D0) THEN
        CI=-1.0D+300 
        SI=0.0D0 
      ELSE IF (X.LE.16.D0) THEN  
        XR=-.25D0*X2 
        CI=EL+DLOG(X)+XR 
        DO K=2,40 
          XR=-. 5D0*XR* (K-1) / (K*K* (2*K-1) ) *X2 
          CI=CI+XR 
          IF (DABS(XR).LT.DABS(CI)*EPS) exit
        enddo
        XR=X 
        SI=X 
        DO K=1,40 
          XR=- . 5D0*XR* (2*K-1) /K/ (4*K*K+4*K+1) *X2 
          SI=SI+XR 
          IF (DABS(XR).LT.DABS(SI)*EPS) RETURN 
        enddo 
      ELSE IF (X.LE.32.0D0) THEN 
        bj=0.d0
        M=INT(47.2+.82*X) 
        XA1=0.0D0 
        XA0=1.0D-100 
        DO K=M,1,-1 
          XA=4.0D0*K*XA0/X-XA1 
          BJ(K)=XA 
          XA1=XA0 
          XA0=XA 
        enddo
        XS=BJ(1) 
        DO K=2,M,2 
          XS=XS+2.0D0*BJ(K+1)
        enddo
        BJ(1)=BJ(1)/XS 
        DO K=1,M 
          BJ(K+1)=BJ(K+1)/XS 
        enddo
        XR=1.0D0 
        XG1=BJ(1) 
        DO K=1,M 
          XR=.25D0*XR*(2*K-1)**2/(K*(2*K+1)**2)*X 
          XG1=XG1+BJ(K+1)*XR 
        enddo
        XR=1.0D0 
        XG2=BJ(1) 
        DO K=1,M 
          XR=.25D0*XR*(2*K-3)**2/(K*(2*K-1)**2)*X 
          XG2=XG2+BJ(K+1)*XR 
        enddo
        XCS=DCOS(X/2.0D0) 
        XSS=DSIN(X/2.0D0) 
        CI=EL+log(X)-X*XSS*XG1+2*XCS*XG2-2*XCS*XCS 
        SI=X*XCS*XG1+2*XSS*XG2-DSIN(X) 
      ELSE
        XR=1.0D0 
        XF=1.0D0 
        DO K=1,9 
          XR=-2.0D0*XR*K*(2*K-1)/X2 
          XF=XF+XR 
        enddo
        XR=1.0D0/X 
        XG=XR 
        DO K=1,8 
          XR=-2.0D0*XR*(2*K+1)*K/X2 
          XG=XG+XR 
        enddo
        CI=XF*DSIN(X)/X-XG*DCOS(X)/X 
        SI=P2-XF*DCOS(X)/X-XG*DSIN(X)/X 
      ENDIF
      deallocate(bj)
      end