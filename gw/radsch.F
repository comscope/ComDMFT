      SUBROUTINE RADSCH(e,z,a,b,l,m,r,dr,p,q,x_q,x_p,v,key1,akap,key,
     &	                p11,q11,p0,q0,key2,nrad,nod1,qq1,work,p1,irel0)
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
c     key  = 1  Old values P and Q are used to start the integration
c	key1 = 1  Exchange is included
c	key2 = 1  Energy derivative finding
c	key2 = 2  Second energy derivative finding
	use manager_mod
	use units_mod
      IMPLICIT none
	integer, intent(in) :: nrad,key,key1,key2,irel0
	real*8, intent(in) :: p1(0:nrad)
	integer :: l,m,nod1,l1,m_rk,mn,i,i0,i00,k,ir
      real*8 :: v(0:nrad),P(0:nrad),Q(0:nrad),R(0:nrad),DR(0:nrad),
     &	      x_p(0:nrad),x_q(0:nrad),q0(0:nrad),p0(0:nrad),
     &          work(0:nrad),arg(10),pp(5),qq(5),e,z,a,b,akap,c2,z2,
     &          dfk2,r3,a3,
     &          b3,v3,xp3,yp3,r1,r2,a1,a2,b1,b2,v1,v2,xp1,xp2,yp1,yp2,
     &          ak1,am1,ak2,am2,ak3,am3,ak4,am4,ac2,ev3,alf,xq3,yq3,
     &          adr1,adr2,adr3,ev1,ev2,xq1,xq2,yq1,yq2,p11,q11,t0,vvv,
     &          ppp,ql,p02,q02,z2c2,alfa,a0,b0,adr,ev,fll,p_p,q_p,qqq,
     &          p_c,q_c,qq1,dqdall
      ppp=0.0d0
      qqq=0.0d0
      xp1=0.0d0
      xp2=0.0d0
      xp3=0.0d0
      yp1=0.0d0
      yp2=0.0d0
      yp3=0.0d0
      L1=L+1
	c2=clight**2
	z2=z+z
	dfk2=dfloat(key2)
	m_rk=5  !! number of Runge-Kutta steps
	mn=6
	do i=1,mn
	  arg(i)=dfloat(i-1)
	enddo
C ****** NON-RELATIVISTIC VERSION (LOGARITHMIC SCALE) ******
      IF(irel0.EQ.0) THEN
C   ASSYMPTOTIC BEHAVIOR ::
        P(0)=0.D0
        Q(0)=0.D0
        I0=1
	  if(key.eq.0) then
   1      P(I0)=R(I0)**L1*(1.D0-z2/(2.D0*L1)*R(I0))
          IF(P(I0).EQ.0.D0)THEN
            Q(I0)=0.D0
            I0=I0+1
            GOTO1
          ENDIF
          Q(I0)=-z2*R(I0)**L1/(L1*2.D0)
	    if(key2.gt.0) then
	      p(i0)=0.d0
	      q(i0)=R(I0)**(L1+1)*z2*z2/(2.D0*L1)/(2*l+3)
	    endif
          IF(Z.LE.1.D-05)Q(I0)=-e*R(I0)**(L1+1)/(2.D0*L1+1)
c         IF(Z.LE.1.D-05)Q(I0)=-(vnuc+EN)*R(I0)**(L1+1)/(2.D0*L1+1)
        else
	    p(i0)=p11
	    q(i0)=q11
        endif
C   4-TH ORDER RUNGE-KUTTA LOOP ::
        R3=R(I0)
        A3=DR(I0)*a
        B3=A3/R3*L1
        V3=A3*(e+v(i0)-z2/r3)
	  if(key1.eq.1) xp3=x_p(i0)*a3
	  if(key2.ge.1) yp3=dfk2*p0(i0)*a3
        DO I=I0,m_rk-1
          R1=R3
          R2=b*(EXP(a*(I+0.5D0))-1.D0) !!! R(I+1/2)
          R3=R(I+1)
          A1=A3
          A2=(R2+B)*a
          A3=DR(I+1)*a
          B1=B3
          B2=A2/R2*L1
          B3=A3/R3*L1
          V1=V3
	    i00=i+1-mn/2
	    i00=max(i00,0)
	    i00=min(i00,nrad-mn+1)
	    t0=dfloat(i-i00)+0.5d0
	    call fint(v(i00),arg,mn,t0,vvv)
          IF(I.NE.nrad-1) V2=A2*
     &		(e+z2/r2-(9.D0*(v(I)+v(I+1))-v(I-1)-v(I+2))*6.25D-2)
          IF(I.EQ.nrad-1) V2=A2*(E+z2/r2-(v(I)+v(I+1))*0.5D0)
	    v2=a2*(e+z2/r2-vvv)
          V3=A3*(E+z2/r3-v(I+1))
	    if(key1.eq.1) then
	      xp1=xp3
            IF(I.NE.nrad-1) xp2=A2*
     &		(9.D0*(x_p(I)+x_p(I+1))-x_p(I-1)-x_p(I+2))*6.25D-2
            IF(I.EQ.nrad-1) xp2=A2*(x_p(I)+x_p(I+1))*0.5D0
	      xp3=x_p(i+1)*a3
	    endif
	    if(key2.ge.1) then
	      yp1=yp3
	    call fint(p0(i00),arg,mn,t0,ppp)
            IF(I.NE.nrad-1) yp2=A2*dfk2*
     &		(9.D0*(p0(I)+p0(I+1))-p0(I-1)-p0(I+2))*6.25D-2
            IF(I.EQ.nrad-1) yp2=dfk2*A2*(p0(I)+p0(I+1))*0.5D0
	      yp2=a2*dfk2*ppp
	      yp3=dfk2*p0(i+1)*a3
	    endif
C   **************
          aK1= A1*Q(I)+B1*P(I)
          aM1=-V1*P(I)-B1*Q(I)
	    if(key1.eq.1) am1=am1+xp1
	    if(key2.ge.1) am1=am1-yp1
C   **************
          aK2= A2*(Q(I)+aM1*0.5D0)+B2*(P(I)+aK1*0.5D0)
          aM2=-V2*(P(I)+aK1*0.5D0)-B2*(Q(I)+aM1*0.5D0)
	    if(key1.eq.1) am2=am2+xp2
	    if(key2.ge.1) am2=am2-yp2
C   **************
          aK3= A2*(Q(I)+aM2*0.5D0)+B2*(P(I)+aK2*0.5D0)
          aM3=-V2*(P(I)+aK2*0.5D0)-B2*(Q(I)+aM2*0.5D0)
	    if(key1.eq.1) am3=am3+xp2
	    if(key2.ge.1) am3=am3-yp2
C   **************
          aK4= A3*(Q(I)+aM3)+B3*(P(I)+aK3)
          aM4=-V3*(P(I)+aK3)-B3*(Q(I)+aM3)
	    if(key1.eq.1) am4=am4+xp3
	    if(key2.ge.1) am4=am4-yp3
C   **************
          P(I+1)=P(I)+(aK1+aK4+2.D0*(aK2+aK3))/6.D0
          Q(I+1)=Q(I)+(aM1+aM4+2.D0*(aM2+aM3))/6.D0
        enddo
C ****** RELATIVISTIC VERSION ******
      ELSEIF(irel0.EQ.1)THEN
C   ASSIMPTOTIC BEHAVIOR ::
        ac2=akap/c2
        QL=SQRT(akap+1.d0-z2**2/C2)
        P(0)=0.D0
        Q(0)=0.D0
        I0=1
	  if(key.eq.0) then
  21      P(I0)=R(I0)**QL
          IF(P(I0).EQ.0.D0)THEN
            Q(I0)=0.D0
            I0=I0+1
            GOTO21
          ENDIF
          Q(I0)=P(I0)*(QL-1)/(Z2/C2)
        else
	    p(i0)=p11
	    q(i0)=q11
	  endif
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
	  if(key2.ge.1) then
	    yp3=dfk2*(1.d0+ac2/(alf*r3)**2)*p0(i0)*adr3
	    yq3=dfk2*q0(i0)*adr3/c2
	  endif
        DO I=I0,m_rk-1
          R1=R3
          R2=B*(EXP(a*(I+0.5D0))-1.D0)   !!! R(I+1/2)
          R3=R(I+1)
	    adr1=adr3
	    adr2=(r2+b)*a
	    adr3=dr(i+1)*a
          B1=B3
          B2=adr2/R2
          B3=adr3/R3
          EV1=EV3
          IF(I.NE.nrad-1) eV2=e+z2/r2
     &		-(9.D0*(v(I)+v(I+1))-v(I-1)-v(I+2))*6.25D-2
          IF(I.EQ.nrad-1) eV2=E+z2/r2-(v(I)+v(I+1))*0.5D0
          EV3=E+z2/r3-v(I+1)
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
	      xp3=x_p(i+1)*adr3
	      xq1=xq3
            IF(I.NE.nrad-1) xq2=adr2/c2*
     &		(9.D0*(x_q(I)+x_q(I+1))-x_q(I-1)-x_q(I+2))*6.25D-2
            IF(I.EQ.nrad-1) xq2=adr2/c2*(x_q(I)+x_q(I+1))*0.5D0
	      xq3=x_q(i+1)*adr3/c2
	    endif
	    if(key2.ge.1) then
	      yp1=yp3
            IF(I.NE.nrad-1) p02=dfk2*
     &		(9.D0*(p0(I)+p0(I+1))-p0(I-1)-p0(I+2))*6.25D-2
            IF(I.EQ.nrad-1) p02=dfk2*(p0(I)+p0(I+1))*0.5D0
	      yp2=(1.d0+ac2/(alf*r2)**2)*p02*adr2
	      yp3=dfk2*(1.d0+ac2/(alf*r3)**2)*p0(i+1)*adr3
	      yq1=yq3
            IF(I.NE.nrad-1) q02=
     &		  (9.D0*(q0(I)+q0(I+1))-q0(I-1)-q0(I+2))*6.25D-2
            IF(I.EQ.nrad-1) q02=(q0(I)+q0(I+1))*0.5D0
	      yq2=dfk2*q02*adr2/c2
	      yq3=dfk2*q0(i+1)*adr3/c2
	    endif
C   **************
          aK1= A1*Q(I)+B1*P(I)
          aM1=-V1*P(I)-B1*Q(I)
	    if(key1.eq.1) then
	      ak1=ak1-xq1
	      am1=am1+xp1
	    endif
	    if(key2.ge.1) then
	      ak1=ak1+yq1
	      am1=am1-yp1
	    endif
C   **************
          aK2= A2*(Q(I)+aM1*0.5D0)+B2*(P(I)+aK1*0.5D0)
          aM2=-V2*(P(I)+aK1*0.5D0)-B2*(Q(I)+aM1*0.5D0)
	    if(key1.eq.1) then
	      ak2=ak2-xq2
	      am2=am2+xp2
	    endif
	    if(key2.ge.1) then
	      ak2=ak2+yq2
	      am2=am2-yp2
	    endif
C   **************
          aK3= A2*(Q(I)+aM2*0.5D0)+B2*(P(I)+aK2*0.5D0)
          aM3=-V2*(P(I)+aK2*0.5D0)-B2*(Q(I)+aM2*0.5D0)
	    if(key1.eq.1) then
	      ak3=ak3-xq2
	      am3=am3+xp2
	    endif
	    if(key2.ge.1) then
	      ak3=ak3+yq2
	      am3=am3-yp2
	    endif
C   **************
          aK4= A3*(Q(I)+aM3)+B3*(P(I)+aK3)
          aM4=-V3*(P(I)+aK3)-B3*(Q(I)+aM3)
	    if(key1.eq.1) then
	      ak4=ak4-xq3
	      am4=am4+xp3
	    endif
	    if(key2.ge.1) then
	      ak4=ak4+yq3
	      am4=am4-yp3
	    endif
C   **************
          P(I+1)=P(I)+(aK1+aK4+2.D0*(aK2+aK3))/6.D0
          Q(I+1)=Q(I)+(aM1+aM4+2.D0*(aM2+aM3))/6.D0
        enddo
      ELSEIF(irel0.EQ.2)THEN
C   ASSIMPTOTIC BEHAVIOR ::
	  z2c2=z2*z2/c2
	  alfa=sqrt(akap*akap-z2c2)
        if(akap.lt.0.d0) then
          a0=1.d0
          if(z.gt.0.5d0) then
            b0=-z2/(alfa-akap)
          else
            b0=0.d0
          endif
        else
          if(z.gt.0.5d0) then
            b0=1.d0
            a0=z2/c2/(alfa+akap)
          else
            a0=0.d0
            b0=1.d0
          endif
        endif
        P(0)=0.D0
        Q(0)=0.D0
        I0=1
	  if(key.eq.0) then
  31      P(I0)=a0*R(I0)**alfa
c          IF(P(I0).EQ.0.D0)THEN
c            Q(I0)=0.D0
c            I0=I0+1
c            GOTO 31
c          ENDIF
          Q(I0)=b0*r(i0)**alfa
        else
	    p(i0)=p11
	    q(i0)=q11
	  endif
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
	  if(key2.ge.1) then
	    yp3=dfk2*p0(i0)*adr3
	    yq3=dfk2*q0(i0)*adr3/c2
	  endif
        DO I=I0,m_rk-1
          R1=R3
          R2=B*(EXP(a*(I+0.5D0))-1.D0)   !!! R(I+1/2)
          R3=R(I+1)
	    adr1=adr3
	    adr2=(r2+b)*a
	    adr3=dr(i+1)*a
          B1=B3
          B2=-akap*adr2/R2
          B3=-akap*adr3/R3
          EV1=EV3
          IF(I.NE.nrad-1) eV2=e+z2/r2
     &		-(9.D0*(v(I)+v(I+1))-v(I-1)-v(I+2))*6.25D-2
          IF(I.EQ.nrad-1) eV2=E+z2/r2-(v(I)+v(I+1))*0.5D0
          EV3=E+z2/r3-v(I+1)
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
	      xp3=x_p(i+1)*adr3
	      xq1=xq3
            IF(I.NE.nrad-1) xq2=adr2/c2*
     &		(9.D0*(x_q(I)+x_q(I+1))-x_q(I-1)-x_q(I+2))*6.25D-2
            IF(I.EQ.nrad-1) xq2=adr2/c2*(x_q(I)+x_q(I+1))*0.5D0
	      xq3=x_q(i+1)*adr3/c2
	    endif
	    if(key2.ge.1) then
	      yp1=yp3
            IF(I.NE.nrad-1) p02=
     &		(9.D0*(p0(I)+p0(I+1))-p0(I-1)-p0(I+2))*6.25D-2
            IF(I.EQ.nrad-1) p02=(p0(I)+p0(I+1))*0.5D0
	      yp2=dfk2*p02*adr2
	      yp3=dfk2*p0(i+1)*adr3
	      yq1=yq3
            IF(I.NE.nrad-1) q02=
     &		  (9.D0*(q0(I)+q0(I+1))-q0(I-1)-q0(I+2))*6.25D-2
            IF(I.EQ.nrad-1) q02=(q0(I)+q0(I+1))*0.5D0
	      yq2=dfk2*q02*adr2/c2
	      yq3=dfk2*q0(i+1)*adr3/c2
	    endif
C   **************
          aK1= A1*Q(I)+B1*P(I)
          aM1=-V1*P(I)-B1*Q(I)
	    if(key1.eq.1) then
	      ak1=ak1-xq1
	      am1=am1+xp1
	    endif
	    if(key2.ge.1) then
	      ak1=ak1+yq1
	      am1=am1-yp1
	    endif
C   **************
          aK2= A2*(Q(I)+aM1*0.5D0)+B2*(P(I)+aK1*0.5D0)
          aM2=-V2*(P(I)+aK1*0.5D0)-B2*(Q(I)+aM1*0.5D0)
	    if(key1.eq.1) then
	      ak2=ak2-xq2
	      am2=am2+xp2
	    endif
	    if(key2.ge.1) then
	      ak2=ak2+yq2
	      am2=am2-yp2
	    endif
C   **************
          aK3= A2*(Q(I)+aM2*0.5D0)+B2*(P(I)+aK2*0.5D0)
          aM3=-V2*(P(I)+aK2*0.5D0)-B2*(Q(I)+aM2*0.5D0)
	    if(key1.eq.1) then
	      ak3=ak3-xq2
	      am3=am3+xp2
	    endif
	    if(key2.ge.1) then
	      ak3=ak3+yq2
	      am3=am3-yp2
	    endif
C   **************
          aK4= A3*(Q(I)+aM3)+B3*(P(I)+aK3)
          aM4=-V3*(P(I)+aK3)-B3*(Q(I)+aM3)
	    if(key1.eq.1) then
	      ak4=ak4-xq3
	      am4=am4+xp3
	    endif
	    if(key2.ge.1) then
	      ak4=ak4+yq3
	      am4=am4-yp3
	    endif
C   **************
          P(I+1)=P(I)+(aK1+aK4+2.D0*(aK2+aK3))/6.D0
          Q(I+1)=Q(I)+(aM1+aM4+2.D0*(aM2+aM3))/6.D0
        enddo
      ENDIF

	
	do ir=1,5
	  i=m_rk+ir-5
	  adr=a*dr(i)
c ------ See V-19 for details -------------------------------------------
	  if(irel0.eq.0) then
		pp(ir)=l1*p(i)/r(i)+q(i)
	    qq(ir)=(v(i)-z2/r(i)-e)*p(i)-l1*q(i)/r(i)
		if(key2.ge.1) qq(ir)=qq(ir)-dfk2*p0(i)
	  else if(irel0.eq.1) then
          ev=e+z2/r(i)-v(i)
	    alf=1.d0+ev/c2
		pp(ir)=p(i)/r(i)+alf*q(i)
	    if(key1.eq.1) pp(ir)=pp(ir)-x_q(i)/c2
	    fll=akap/r(i)**2/alf
	    qq(ir)=(fll-ev)*p(i)-q(i)/r(i)
		if(key2.ge.1) then
	      pp(ir)=pp(ir)+dfk2*q0(i)/c2
		  qq(ir)=qq(ir)-dfk2*(1.d0+fll/alf/c2)*p0(i)
		  if(key2.eq.2) qq(ir)=qq(ir)+2.d0*fll/alf**2/c2**2*p1(i)
          endif
	  else if(irel0.eq.2) then
          ev=e+z2/r(i)-v(i)
	    alf=1.d0+ev/c2
		pp(ir)=-akap*p(i)/r(i)+alf*q(i)
	    if(key1.eq.1) pp(ir)=pp(ir)-x_q(i)/c2
	    qq(ir)=-ev*p(i)+akap*q(i)/r(i)
		if(key2.ge.1) then
	      pp(ir)=pp(ir)+dfk2*q0(i)/c2
		  qq(ir)=qq(ir)-dfk2*p0(i)
          endif
	  endif
	  pp(ir)=adr*pp(ir)
	  if(key1.eq.1) qq(ir)=qq(ir)+x_p(i)
	  qq(ir)=adr*qq(ir)
	enddo

	do i=m_rk,m-1
	  p_p=p(i)+(1901.d0*pp(5)-2774.d0*pp(4)+2616.d0*pp(3)
     &	       -1274.d0*pp(2)+251.d0*pp(1))/720.d0
	  q_p=q(i)+(1901.d0*qq(5)-2774.d0*qq(4)+2616.d0*qq(3)
     &	       -1274.d0*qq(2)+251.d0*qq(1))/720.d0
	  adr=a*dr(i+1)
c ------ See V-19 for details -------------------------------------------
	  if(irel0.eq.0) then
		  ppp=l1*p_p/r(i+1)+q_p
	    qqq=(v(i+1)-z2/r(i+1)-e)*p_p-l1*q_p/r(i+1)
		  if(key2.ge.1) qqq=qqq-dfk2*p0(i+1)
	  else if(irel0.eq.1) then
          ev=e+z2/r(i+1)-v(i+1)
	    alf=1.d0+ev/c2
		  ppp=p_p/r(i+1)+q_p*alf
	    if(key1.eq.1) ppp=ppp-x_q(i+1)/c2
	    fll=akap/r(i+1)**2/alf
	    qqq=(fll-ev)*p_p-q_p/r(i+1)
		  if(key2.ge.1) then
	      ppp=ppp+dfk2*q0(i+1)/c2
		    qqq=qqq-dfk2*(1.d0+fll/alf/c2)*p0(i+1)
		    if(key2.eq.2) qqq=qqq+2.d0*fll/alf**2/c2**2*p1(i+1)
          endif
	  else if(irel0.eq.2) then
          ev=e+z2/r(i+1)-v(i+1)
	    alf=1.d0+ev/c2
		  ppp=-akap*p_p/r(i+1)+alf*q_p
	    if(key1.eq.1) ppp=ppp-x_q(i+1)/c2
	    qqq=-ev*p_p+akap*q_p/r(i+1)
		  if(key2.ge.1) then
	      ppp=ppp+dfk2*q0(i+1)/c2
		    qqq=qqq-dfk2*p0(i+1)
          endif
	  endif
	  ppp=adr*ppp
	  if(key1.eq.1) qqq=qqq+x_p(i+1)
	  qqq=adr*qqq
	  p_c=p(i)+(251.d0*ppp+646.d0*pp(5)-264.d0*pp(4)
     &	       +106.d0*pp(3)-19.d0*pp(2))/720.d0
	  p(i+1)=(475.d0*p_c+27.d0*p_p)/502.d0
	  q_c=q(i)+(251.d0*qqq+646.d0*qq(5)-264.d0*qq(4)
     &	       +106.d0*qq(3)-19.d0*qq(2))/720.d0
	  q(i+1)=(475.d0*q_c+27.d0*q_p)/502.d0
	  do k=1,4
	    pp(k)=pp(k+1)
	    qq(k)=qq(k+1)
	  enddo
	  adr=a*dr(i+1)
c ------ See V-19 for details -------------------------------------------
	  if(irel0.eq.0) then
		  pp(5)=l1*p(i+1)/r(i+1)+q(i+1)
	    qq(5)=(v(i+1)-z2/r(i+1)-e)*p(i+1)-l1*q(i+1)/r(i+1)
		  if(key2.ge.1) qq(5)=qq(5)-dfk2*p0(i+1)
	  else if(irel0.eq.1) then
          ev=e+z2/r(i+1)-v(i+1)
	    alf=1.d0+ev/c2
		  pp(5)=p(i+1)/r(i+1)+q(i+1)*alf
	    if(key1.eq.1) pp(5)=pp(5)-x_q(i+1)/c2
	    fll=akap/r(i+1)**2/alf
	    qq(5)=(fll-ev)*p(i+1)-q(i+1)/r(i+1)
		  if(key2.ge.1) then
	      pp(5)=pp(5)+dfk2*q0(i+1)/c2
		    qq(5)=qq(5)-dfk2*(1.d0+fll/alf/c2)*p0(i+1)
		    if(key2.eq.2) qq(5)=qq(5)+2.d0*fll/alf**2/c2**2*p1(i+1)
          endif
	  else if(irel0.eq.2) then
          ev=e+z2/r(i+1)-v(i+1)
	    alf=1.d0+ev/c2
		  pp(5)=-akap*p(i+1)/r(i+1)+alf*q(i+1)
	    if(key1.eq.1) pp(5)=pp(5)-x_q(i+1)/c2
	    qq(5)=-ev*p(i+1)+akap*q(i+1)/r(i+1)
		  if(key2.ge.1) then
	      pp(5)=pp(5)+dfk2*q0(i+1)/c2
		    qq(5)=qq(5)-dfk2*p0(i+1)
          endif
	  endif
	  pp(5)=adr*pp(5)
	  if(key1.eq.1) qq(5)=qq(5)+x_p(i+1)
	  qq(5)=adr*qq(5)
	enddo

	nod1=0
	do ir=0,m-1
	  if(p(ir)*p(ir+1).lt.0.d0) nod1=nod1+1
	enddo

	do ir=0,m
	  work(ir)=p(ir)**2*dr(ir)
      enddo
	if(irel0.ge.1) then
	  do ir=0,m
	    work(ir)=work(ir)+q(ir)**2*dr(ir)/c2
        enddo
	endif
      qq1=dqdall(a,work(0),m)
      END
