        SUBROUTINE RADSCH_b(e,z,a,b,l,m,r,dr,p,q,x_q,x_p,v,key1,akap,
     &	                    p1,q1,p0,q0,key2,nrad,nod2,q2,work,irel0)
***************************************************************
*  The solution of the radial Schroedinger equation           *
*                                                             *
*   / P'(r)=Q(r)+(l+1)*P(r)/r                                 *
*   \ Q'(r)=(V(r)-E)*P(r)-(l+1)*Q(r)/r  + Xp                  *
*     FI(r)=P(r)/r and FID(r)=FI'(r)                          *
*                                                             *
*  or scalar-relativistic Dirac equation                      *
*                                                             *
*   / P'(r)=A(r)*Q(r)+P(r)/r-Xq/c^2                           *
*   \ Q'(r)=(V(r)-E+C(r))*P(r)-Q(r)/r + Xp                    *
*  where                                                      *
*     A(r)=1+(E-V(r))/c^2                                     *
*     C(r)=l*(l+1)/r^2/A(r)                                   *
*     FI(r)=P(r)/r and FID(r)=FI'(r)                          *
*  and c is the light velocity                                *
***************************************************************
c	key1 = 1  Exchange is included
c	key2 = 1  Energy derivative finding
	use manager_mod
	use units_mod
      IMPLICIT none
      integer, intent(in) :: irel0,nrad
	integer :: l,m,key1,key2,nod2,l1,i0,i,kc,ir,ii,k
      real*8 :: v(0:nrad),P(0:nrad),Q(0:nrad),R(0:nrad),DR(0:nrad),
     &	        x_p(0:nrad),x_q(0:nrad),q0(0:nrad),p0(0:nrad),
     &          work(0:nrad),pp(5),qq(5),e,z,a,b,akap,c2,z2,dfk2,r3,a3,
     &          b3,v3,xp3,yp3,r1,r2,a1,a2,b1,b2,v1,v2,xp1,xp2,yp1,yp2,
     &          ak1,am1,ak2,am2,ak3,am3,ak4,am4,ac2,ev3,alf,xq3,yq3,
     &          adr1,adr2,adr3,ev1,ev2,xq1,xq2,yq1,yq2,q2,p1,q1,p02,q02,
     &          adr,ev,fll,p_p,q_p,ppp,qqq,p_c,q_c,dqdall
	c2=clight**2
      ppp=0.0d0
      qqq=0.0d0
      xp1=0.0d0
      xp2=0.0d0
      xp3=0.0d0
      yp1=0.0d0
      yp2=0.0d0
      yp3=0.0d0
      L1=L+1
	z2=z+z
	dfk2=dfloat(key2)
C ****** NON-RELATIVISTIC VERSION (LOGARITHMIC SCALE) ******
      IF(irel0.EQ.0) THEN
C   ASYMPTOTIC BEHAVIOR ::
        I0=nrad
	  p(i0)=p1
	  q(i0)=q1
C   4-TH ORDER RUNGE-KUTTA LOOP ::
        R3=R(I0)
        A3=DR(I0)*a
        B3=A3/R3*L1
        V3=A3*(e+v(i0)-z2/r3)
	  if(key1.eq.1) xp3=x_p(i0)*a3
	  if(key2.eq.1) yp3=p0(i0)*a3
        DO I=nrad-1,nrad-4,-1
          R1=R3
          R2=b*(EXP(a*(I+0.5D0))-1.D0) !!! R(I+1/2)
          R3=R(I)
          A1=A3
          A2=(R2+B)*a
          A3=DR(I)*a
          B1=B3
          B2=A2/R2*L1
          B3=A3/R3*L1
          V1=V3
          IF(I.NE.nrad-1) V2=A2*
     &		(e+z2/r2-(9.D0*(v(I)+v(I+1))-v(I-1)-v(I+2))*6.25D-2)
          IF(I.EQ.nrad-1) V2=A2*(E+z2/r2-(v(I)+v(I+1))*0.5D0)
          V3=A3*(E+z2/r3-v(I))
	    if(key1.eq.1) then
	      xp1=xp3
            IF(I.NE.nrad-1) xp2=A2*
     &		(9.D0*(x_p(I)+x_p(I+1))-x_p(I-1)-x_p(I+2))*6.25D-2
            IF(I.EQ.nrad-1) xp2=A2*(x_p(I)+x_p(I+1))*0.5D0
	      xp3=x_p(i)*a3
	    endif
	    if(key2.eq.1) then
	      yp1=yp3
            IF(I.NE.nrad-1) yp2=A2*
     &		(9.D0*(p0(I)+p0(I+1))-p0(I-1)-p0(I+2))*6.25D-2
            IF(I.EQ.nrad-1) yp2=A2*(p0(I)+p0(I+1))*0.5D0
	      yp3=p0(i)*a3
	    endif
C   **************
          aK1= A1*Q(I+1)+B1*P(I+1)
          aM1=-V1*P(I+1)-B1*Q(I+1)
	    if(key1.eq.1) am1=am1+xp1
	    if(key2.eq.1) am1=am1-yp1
C   **************
          aK2= A2*(Q(I+1)-aM1*0.5D0)+B2*(P(I+1)-aK1*0.5D0)
          aM2=-V2*(P(I+1)-aK1*0.5D0)-B2*(Q(I+1)-aM1*0.5D0)
	    if(key1.eq.1) am2=am2+xp2
	    if(key2.eq.1) am2=am2-yp2
C   **************
          aK3= A2*(Q(I+1)-aM2*0.5D0)+B2*(P(I+1)-aK2*0.5D0)
          aM3=-V2*(P(I+1)-aK2*0.5D0)-B2*(Q(I+1)-aM2*0.5D0)
	    if(key1.eq.1) am3=am3+xp2
	    if(key2.eq.1) am3=am3-yp2
C   **************
          aK4= A3*(Q(I+1)-aM3)+B3*(P(I+1)-aK3)
          aM4=-V3*(P(I+1)-aK3)-B3*(Q(I+1)-aM3)
	    if(key1.eq.1) am4=am4+xp3
	    if(key2.eq.1) am4=am4-yp3
C   **************
          P(I)=P(I+1)-(aK1+aK4+2.D0*(aK2+aK3))/6.D0
          Q(I)=Q(I+1)-(aM1+aM4+2.D0*(aM2+aM3))/6.D0
        enddo
C ****** RELATIVISTIC VERSION ******
      ELSEIF(irel0.EQ.1)THEN
C   ASyMPTOTIC BEHAVIOR ::
        ac2=akap/c2
        I0=nrad
	  p(i0)=p1
	  q(i0)=q1
C   4-TH ORDER RUNGE-KUTTA LOOP ::
        R3=R(I0)
	  adr3=a*dr(i0)
        B3=adr3/R3
        EV3=E+z2/R3-v(i0)
	  alf=1.d0+ev3/c2
        A3=adr3*alf
        V3=adr3*(EV3-akap/R3/R3/alf)
	  if(key1.eq.1) then
	    xp3=x_p(i0)*adr3
	    xq3=x_q(i0)*adr3/c2
	  endif
	  if(key2.eq.1) then
	    yp3=(1.d0+ac2/(alf*r3)**2)*p0(i0)*adr3
	    yq3=q0(i0)*adr3/c2
	  endif
        DO I=nrad-1,nrad-4,-1
          R1=R3
          R2=B*(EXP(a*(I+0.5D0))-1.D0)   !!! R(I+1/2)
          R3=R(I)
	    adr1=adr3
	    adr2=(r2+b)*a
	    adr3=dr(i)*a
          B1=B3
          B2=adr2/R2
          B3=adr3/R3
          EV1=EV3
          IF(I.NE.nrad-1) eV2=e+z2/r2
     &		-(9.D0*(v(I)+v(I+1))-v(I-1)-v(I+2))*6.25D-2
          IF(I.EQ.nrad-1) eV2=E+z2/r2-(v(I)+v(I+1))*0.5D0
          EV3=E+z2/r3-v(I)
          A1=A3
          A2=adr2*(1.D0+EV2/C2)
          A3=adr3*(1.D0+EV3/C2)
          V1=V3
          V2=adr2*(EV2-akap/R2/R2/(1.D0+EV2/C2))
          V3=adr3*(EV3-akap/R3/R3/(1.D0+EV3/C2))
	    if(key1.eq.1) then
	      xp1=xp3
            IF(I.NE.nrad-1) xp2=adr2*
     &		(9.D0*(x_p(I)+x_p(I+1))-x_p(I-1)-x_p(I+2))*6.25D-2
            IF(I.EQ.nrad-1) xp2=adr2*(x_p(I)+x_p(I+1))*0.5D0
	      xp3=x_p(i)*adr3
	      xq1=xq3
            IF(I.NE.nrad-1) xq2=adr2/c2*
     &		(9.D0*(x_q(I)+x_q(I+1))-x_q(I-1)-x_q(I+2))*6.25D-2
            IF(I.EQ.nrad-1) xq2=adr2/c2*(x_q(I)+x_q(I+1))*0.5D0
	      xq3=x_q(i)*adr3/c2
	    endif
	    if(key2.eq.1) then
	      yp1=yp3
            IF(I.NE.nrad-1) p02=
     &		(9.D0*(p0(I)+p0(I+1))-p0(I-1)-p0(I+2))*6.25D-2
            IF(I.EQ.nrad-1) p02=(p0(I)+p0(I+1))*0.5D0
	      yp2=(1.d0+ac2/(alf*r2)**2)*p02*adr2
	      yp3=(1.d0+ac2/(alf*r3)**2)*p0(i)*adr3
	      yq1=yq3
            IF(I.NE.nrad-1) q02=
     &		  (9.D0*(q0(I)+q0(I+1))-q0(I-1)-q0(I+2))*6.25D-2
            IF(I.EQ.nrad-1) q02=(q0(I)+q0(I+1))*0.5D0
	      yq2=q02*adr2/c2
	      yq3=q0(i)*adr3/c2
	    endif
C   **************
          aK1= A1*Q(I+1)+B1*P(I+1)
          aM1=-V1*P(I+1)-B1*Q(I+1)
	    if(key1.eq.1) then
	      ak1=ak1-xq1
	      am1=am1+xp1
	    endif
	    if(key2.eq.1) then
	      ak1=ak1+yq1
	      am1=am1-yp1
	    endif
C   **************
          aK2= A2*(Q(I+1)-aM1*0.5D0)+B2*(P(I+1)-aK1*0.5D0)
          aM2=-V2*(P(I+1)-aK1*0.5D0)-B2*(Q(I+1)-aM1*0.5D0)
	    if(key1.eq.1) then
	      ak2=ak2-xq2
	      am2=am2+xp2
	    endif
	    if(key2.eq.1) then
	      ak2=ak2+yq2
	      am2=am2-yp2
	    endif
C   **************
          aK3= A2*(Q(I+1)-aM2*0.5D0)+B2*(P(I+1)-aK2*0.5D0)
          aM3=-V2*(P(I+1)-aK2*0.5D0)-B2*(Q(I+1)-aM2*0.5D0)
	    if(key1.eq.1) then
	      ak3=ak3-xq2
	      am3=am3+xp2
	    endif
	    if(key2.eq.1) then
	      ak3=ak3+yq2
	      am3=am3-yp2
	    endif
C   **************
          aK4= A3*(Q(I+1)-aM3)+B3*(P(I+1)-aK3)
          aM4=-V3*(P(I+1)-aK3)-B3*(Q(I+1)-aM3)
	    if(key1.eq.1) then
	      ak4=ak4-xq3
	      am4=am4+xp3
	    endif
	    if(key2.eq.1) then
	      ak4=ak4+yq3
	      am4=am4-yp3
	    endif
C   **************
          P(I)=P(I+1)-(aK1+aK4+2.D0*(aK2+aK3))/6.D0
          Q(I)=Q(I+1)-(aM1+aM4+2.D0*(aM2+aM3))/6.D0
        enddo
      ELSEIF(irel0.EQ.2)THEN
C   ASyMPTOTIC BEHAVIOR ::
        I0=nrad
	  p(i0)=p1
	  q(i0)=q1
C   4-TH ORDER RUNGE-KUTTA LOOP ::
        R3=R(I0)
	  adr3=a*dr(i0)
        B3=-akap*adr3/R3
        EV3=E+z2/R3-v(i0)
	  alf=1.d0+ev3/c2
        A3=adr3*alf
        V3=adr3*EV3
	  if(key1.eq.1) then
	    xp3=x_p(i0)*adr3
	    xq3=x_q(i0)*adr3/c2
	  endif
	  if(key2.eq.1) then
	    yp3=p0(i0)*adr3
	    yq3=q0(i0)*adr3/c2
	  endif
        DO I=nrad-1,nrad-4,-1
          R1=R3
          R2=B*(EXP(a*(I+0.5D0))-1.D0)   !!! R(I+1/2)
          R3=R(I)
	    adr1=adr3
	    adr2=(r2+b)*a
	    adr3=dr(i)*a
          B1=B3
          B2=-akap*adr2/R2
          B3=-akap*adr3/R3
          EV1=EV3
          IF(I.NE.nrad-1) eV2=e+z2/r2
     &		-(9.D0*(v(I)+v(I+1))-v(I-1)-v(I+2))*6.25D-2
          IF(I.EQ.nrad-1) eV2=E+z2/r2-(v(I)+v(I+1))*0.5D0
          EV3=E+z2/r3-v(I)
          A1=A3
          A2=adr2*(1.D0+EV2/C2)
          A3=adr3*(1.D0+EV3/C2)
          V1=V3
          V2=adr2*EV2
          V3=adr3*EV3
	    if(key1.eq.1) then
	      xp1=xp3
            IF(I.NE.nrad-1) xp2=adr2*
     &		(9.D0*(x_p(I)+x_p(I+1))-x_p(I-1)-x_p(I+2))*6.25D-2
            IF(I.EQ.nrad-1) xp2=adr2*(x_p(I)+x_p(I+1))*0.5D0
	      xp3=x_p(i)*adr3
	      xq1=xq3
            IF(I.NE.nrad-1) xq2=adr2/c2*
     &		(9.D0*(x_q(I)+x_q(I+1))-x_q(I-1)-x_q(I+2))*6.25D-2
            IF(I.EQ.nrad-1) xq2=adr2/c2*(x_q(I)+x_q(I+1))*0.5D0
	      xq3=x_q(i)*adr3/c2
	    endif
	    if(key2.eq.1) then
	      yp1=yp3
            IF(I.NE.nrad-1) p02=
     &		(9.D0*(p0(I)+p0(I+1))-p0(I-1)-p0(I+2))*6.25D-2
            IF(I.EQ.nrad-1) p02=(p0(I)+p0(I+1))*0.5D0
	      yp2=p02*adr2
	      yp3=p0(i)*adr3
	      yq1=yq3
            IF(I.NE.nrad-1) q02=
     &		  (9.D0*(q0(I)+q0(I+1))-q0(I-1)-q0(I+2))*6.25D-2
            IF(I.EQ.nrad-1) q02=(q0(I)+q0(I+1))*0.5D0
	      yq2=q02*adr2/c2
	      yq3=q0(i)*adr3/c2
	    endif
C   **************
          aK1= A1*Q(I+1)+B1*P(I+1)
          aM1=-V1*P(I+1)-B1*Q(I+1)
	    if(key1.eq.1) then
	      ak1=ak1-xq1
	      am1=am1+xp1
	    endif
	    if(key2.eq.1) then
	      ak1=ak1+yq1
	      am1=am1-yp1
	    endif
C   **************
          aK2= A2*(Q(I+1)-aM1*0.5D0)+B2*(P(I+1)-aK1*0.5D0)
          aM2=-V2*(P(I+1)-aK1*0.5D0)-B2*(Q(I+1)-aM1*0.5D0)
	    if(key1.eq.1) then
	      ak2=ak2-xq2
	      am2=am2+xp2
	    endif
	    if(key2.eq.1) then
	      ak2=ak2+yq2
	      am2=am2-yp2
	    endif
C   **************
          aK3= A2*(Q(I+1)-aM2*0.5D0)+B2*(P(I+1)-aK2*0.5D0)
          aM3=-V2*(P(I+1)-aK2*0.5D0)-B2*(Q(I+1)-aM2*0.5D0)
	    if(key1.eq.1) then
	      ak3=ak3-xq2
	      am3=am3+xp2
	    endif
	    if(key2.eq.1) then
	      ak3=ak3+yq2
	      am3=am3-yp2
	    endif
C   **************
          aK4= A3*(Q(I+1)-aM3)+B3*(P(I+1)-aK3)
          aM4=-V3*(P(I+1)-aK3)-B3*(Q(I+1)-aM3)
	    if(key1.eq.1) then
	      ak4=ak4-xq3
	      am4=am4+xp3
	    endif
	    if(key2.eq.1) then
	      ak4=ak4+yq3
	      am4=am4-yp3
	    endif
C   **************
          P(I)=P(I+1)-(aK1+aK4+2.D0*(aK2+aK3))/6.D0
          Q(I)=Q(I+1)-(aM1+aM4+2.D0*(aM2+aM3))/6.D0
        enddo
      ENDIF


	do i=1,5
	  ii=nrad-i+1
	  adr=a*dr(ii)
c ------ See V-19 for details -------------------------------------------
	  if(irel0.eq.0) then
	    pp(i)=l1*p(ii)/r(ii)+q(ii)
	    qq(i)=(v(ii)-z2/r(ii)-e)*p(ii)-l1*q(ii)/r(ii)
		  if(key2.ge.1) qq(i)=qq(i)-dfk2*p0(ii)
	  else if(irel0.eq.1) then
          ev=e+z2/r(ii)-v(ii)
	    alf=1.d0+ev/c2
		  pp(i)=p(ii)/r(ii)+alf*q(ii)
	    if(key1.eq.1) pp(i)=pp(i)-x_q(ii)/c2
	    fll=akap/r(ii)**2/alf
	    qq(i)=(fll-ev)*p(ii)-q(ii)/r(ii)
		  if(key2.ge.1) then
	      pp(i)=pp(i)+dfk2*q0(ii)/c2
		    qq(i)=qq(i)-dfk2*(1.d0+fll/alf/c2)*p0(ii)
c		    if(key2.eq.2) qq(i)=qq(i)+2.d0*fll/alf**2/c2**2*p1(ii)
          endif
	  else if(irel0.eq.2) then
          ev=e+z2/r(ii)-v(ii)
	    alf=1.d0+ev/c2
		  pp(i)=-akap*p(ii)/r(ii)+alf*q(ii)
	    if(key1.eq.1) pp(i)=pp(i)-x_q(ii)/c2
	    qq(i)=-ev*p(ii)+akap*q(ii)/r(ii)
		  if(key2.ge.1) then
	      pp(i)=pp(i)+dfk2*q0(ii)/c2
		    qq(i)=qq(i)-dfk2*p0(ii)
          endif
	  endif
	  pp(i)=adr*pp(i)
	  if(key1.eq.1) qq(i)=qq(i)+x_p(ii)
	  qq(i)=adr*qq(i)
	enddo

	do i=nrad-4,m+1,-1
	  p_p=p(i)-(1901.d0*pp(5)-2774.d0*pp(4)+2616.d0*pp(3)
     &	       -1274.d0*pp(2)+251.d0*pp(1))/720.d0
	  q_p=q(i)-(1901.d0*qq(5)-2774.d0*qq(4)+2616.d0*qq(3)
     &	       -1274.d0*qq(2)+251.d0*qq(1))/720.d0
	  adr=a*dr(i-1)
c ------ See V-19 for details -------------------------------------------
	  if(irel0.eq.0) then
		  ppp=l1*p_p/r(i-1)+q_p
	    qqq=(v(i-1)-z2/r(i-1)-e)*p_p-l1*q_p/r(i-1)
		  if(key2.ge.1) qqq=qqq-dfk2*p0(i-1)
	  else if(irel0.eq.1) then
          ev=e+z2/r(i-1)-v(i-1)
	    alf=1.d0+ev/c2
		  ppp=p_p/r(i-1)+alf*q_p
	    if(key1.eq.1) ppp=ppp-x_q(i-1)/c2
	    fll=akap/r(i-1)**2/alf
	    qqq=(fll-ev)*p_p-q_p/r(i-1)
		  if(key2.ge.1) then
	      ppp=ppp+dfk2*q0(i-1)/c2
		    qqq=qqq-dfk2*(1.d0+fll/alf/c2)*p0(i-1)
c		    if(key2.eq.2) qqq=qqq+2.d0*fll/alf**2/c2**2*p1(i-1)
          endif
	  else if(irel0.eq.2) then
          ev=e+z2/r(i-1)-v(i-1)
	    alf=1.d0+ev/c2
		  ppp=-akap*p_p/r(i-1)+alf*q_p
	    if(key1.eq.1) ppp=ppp-x_q(i-1)/c2
	    qqq=-ev*p_p+akap*q_p/r(i-1)
		  if(key2.ge.1) then
	      ppp=ppp+dfk2*q0(i-1)/c2
		    qqq=qqq-dfk2*p0(i-1)
          endif
	  endif
	  ppp=adr*ppp
	  if(key1.eq.1) qqq=qqq+x_p(i-1)
	  qqq=adr*qqq
	  p_c=p(i)-(251.d0*ppp+646.d0*pp(5)-264.d0*pp(4)
     &	        +106.d0*pp(3)-19.d0*pp(2))/720.d0
	  p(i-1)=(475.d0*p_c+27.d0*p_p)/502.d0
	  q_c=q(i)-(251.d0*qqq+646.d0*qq(5)-264.d0*qq(4)
     &	        +106.d0*qq(3)-19.d0*qq(2))/720.d0
	  q(i-1)=(475.d0*q_c+27.d0*q_p)/502.d0
	  do k=1,4
	    pp(k)=pp(k+1)
	    qq(k)=qq(k+1)
	  enddo
	  adr=a*dr(i-1)
c ------ See V-19 for details -------------------------------------------
	  if(irel0.eq.0) then
		  pp(5)=l1*p(i-1)/r(i-1)+q(i-1)
	    qq(5)=(v(i-1)-z2/r(i-1)-e)*p(i-1)-l1*q(i-1)/r(i-1)
		  if(key2.ge.1) qq(5)=qq(5)-dfk2*p0(i-1)
	  else if(irel0.eq.1) then
          ev=e+z2/r(i-1)-v(i-1)
	    alf=1.d0+ev/c2
		  pp(5)=p(i-1)/r(i-1)+alf*q(i-1)
	    if(key1.eq.1) pp(5)=pp(5)-x_q(i-1)/c2
	    fll=akap/r(i-1)**2/alf
	    qq(5)=(fll-ev)*p(i-1)-q(i-1)/r(i-1)
		  if(key2.ge.1) then
	      pp(5)=pp(5)+dfk2*q0(i-1)/c2
		    qq(5)=qq(5)-dfk2*(1.d0+fll/alf/c2)*p0(i-1)
c		    if(key2.eq.2) qq(5)=qq(5)+2.d0*fll/alf**2/c2**2*p1(i-1)
          endif
	  else if(irel0.eq.2) then
          ev=e+z2/r(i-1)-v(i-1)
	    alf=1.d0+ev/c2
		  pp(5)=-akap*p(i-1)/r(i-1)+alf*q(i-1)
	    if(key1.eq.1) pp(5)=pp(5)-x_q(i-1)/c2
	    qq(5)=-ev*p(i-1)+akap*q(i-1)/r(i-1)
		  if(key2.ge.1) then
	      pp(5)=pp(5)+dfk2*q0(i-1)/c2
		    qq(5)=qq(5)-dfk2*p0(i-1)
          endif
	  endif
	  pp(5)=adr*pp(5)
	  if(key1.eq.1) qq(5)=qq(5)+x_p(i-1)
	  qq(5)=adr*qq(5)
	  if(key1.eq.0.and.mod(i-1,2).eq.0) then
		  if(p(i-1)*pp(5).gt.0.d0) then
	      kc=i-1
            goto 2
	    endif
	  endif
	enddo
	kc=m
2	continue
	m=kc



	nod2=0
	do ir=m,nrad-1
	  if(p(ir)*p(ir+1).lt.0.d0) nod2=nod2+1
	enddo

	do ir=m,nrad
	  work(ir)=p(ir)**2*dr(ir)
      enddo
	if(irel0.ge.1) then
	  do ir=m,nrad
	    work(ir)=work(ir)+q(ir)**2*dr(ir)/c2
        enddo
	endif
      q2=dqdall(a,work(m),nrad-m)
      END
