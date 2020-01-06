

      subroutine rad_eq(eb1,eb2,e,tol,z,l,i,nod,val,slo,v,q,r,nrad,pnew,
     &  qnew,f_ex,fi00,dfi00,kc,nit,de,nre,f_exq,dr,
     &  key1,h,dl,work,acc,node,irel0,clight)
      
c     The Bound solution for given L and NOD
c     key1 = 1  Exchange is included
c     Input:
c     val = P(S)
c     slo = P'(S)
c     nod = number of nodes
c     f_ex =
c     f_exq = 
c     Output:
c     pnew = P_new(r)
c     qnew = Q_new(r)
c     fi =  P(S)/S
c     dfi = (P(S)/S)'
c     dl = log derivative
c     q = <P!P>
      implicit none
      
      integer, intent(in) :: irel0,nrad,nod,key1,l,i
      double precision, intent(in) :: v(0:nrad),val,slo,h,r(0:nrad),
     $  dr(0:nrad),f_exq(0:nrad),f_ex(0:nrad),z,tol,eb1,eb2,
     $  clight
      integer, intent(inout) :: nre
      double precision, intent(inout) :: e,pnew(0:nrad),qnew(0:nrad)
      integer, intent(out) :: node,kc,nit
      double precision, intent(out) :: acc,dl,de,dfi00,fi00,q,
     $  work(0:nrad)

      integer :: nitmax,k,nod2,nod1,ir
      double precision :: c2,e1,e2,akap,
     &  pn,qn,alf,ph_in,qh_in,ratio,q1,q2,pi_in,qi_in,e00,
     &  q001,q01,ee,fff,b,p0,q0,p1,qq1,e0,fac,dp,dfi_big,h12,
     &  emi,ema
      logical bis,tw
      double precision, allocatable :: g1(:),g2(:)

      nitmax=5000
      c2=clight**2
      q001=1.0d15
      de=1.0d0
      nit=0
c     if(l==0.and.nod==2) nit=1
      e1=eb1
      e2=eb2
      e00=e
      b=r(nrad)/(dexp(h*nrad-h)-1.d0)
c     ----- start iterations to find energy ------
      bis=.false.
      tw=.false.
      emi=e
      ema=e  
 1    nit=nit+1
      if(nit.gt.nitmax) return
      if(e.le.e1.or.e.ge.e2) e=.5d0*(e1+e2)
      emi=min(emi,e)
      ema=max(ema,e)
      if(key1==0.or.nre==nrad) then
        call get_nre(e,nre,l,r,v,z,nrad)
        nre=min(nre,nrad)
      endif
      call pn_qn(irel0,e,z,l,i,val,slo,v(nre),r(nre),pn,qn,f_exq(nre),
     &  key1,akap,clight)
      kc=nre*2/3
      if(mod(kc,2).eq.1) kc=kc-1
      if(kc.ge.nre-20) kc=nre-20
c     -------- Homogeneous equation ---------------------------------------
c     ------- Inward integration ----------------------------------
      call RADSCH_b(e,z,h,b,l,kc,r,dr,pnew,qnew,f_exq,f_ex,v,0,
     &  akap,pn,qn,p0,q0,0,nre,nod2,q2,work,irel0,clight)
      ph_in=pnew(kc)
      qh_in=qnew(kc)
c     Outward integration ----------------------------------------------
      call RADSCH(e,z,h,b,l,kc,r,dr,pnew,qnew,f_exq,f_ex,v,0,
     &  akap,0,p1,qq1,p0,q0,0,nre,nod1,q1,work,p0,irel0,clight)
      if(key1.eq.0) then        !! id est exchange is not included
        ratio=pnew(kc)/ph_in
        do ir=0,kc
          pnew(ir)=pnew(ir)/ratio
          qnew(ir)=qnew(ir)/ratio
        enddo
        q1=q1/ratio**2
        q=q1+q2
        node=nod1+nod2
c     call rad_hf_check(pnew,qnew,nre,node,nrad)
      else if(key1.eq.1) then   !! id est exchange is included
c     -------- Inhomogeneous equation --------------------------------------
        allocate(g1(0:nre))
        allocate(g2(0:nre))
c     ------- Inward integration ----------------------------------
        call RADSCH_b(e,z,h,b,l,kc,r,dr,g1,g2,f_exq,f_ex,v,1,
     &    akap,pn,qn,p0,q0,0,nre,nod2,q2,work,irel0,clight)
     	  pi_in=g1(kc)
     	  qi_in=g2(kc)
c     Outward integration ----------------------------------------------
        call RADSCH(e,z,h,b,l,kc,r,dr,g1,g2,f_exq,f_ex,v,1,
     &    akap,0,p1,qq1,p0,q0,0,nre,nod1,q1,work,p0,irel0,clight)
        call inhomo(ph_in,pi_in,qh_in,qi_in,pnew,qnew,g1,g2,
     &    node,q,nre,kc,dr,h,irel0,clight)
        deallocate(g1,g2)
        call rad_hf_check(pnew,qnew,nre,node,nrad)
      endif
      if(key1==0) then
        if(node.ne.nod) then
          if(node.gt.nod) e2=e
          if(node.lt.nod) e1=e
          if(node.gt.nod) e=min(e*1.5d0,e-0.1d0)
          if(node.lt.nod) e=max(e*0.7d0,e+0.1d0)
          goto 1
        endif
      endif
      if(key1.eq.0) then
        de=ph_in*(qnew(kc)-qh_in)/q
      else
        if(.not.bis) then
          if(.not.tw) then
            de=-max(0.1d0,abs(e)/2000.d0)
            if(node<nod) de=-de
            e00=e
            q001=q-1.d0
            tw=.true.
          else
            e0=e
            q01=q-1.d0
            if(q01*q001.lt.0.d0) then
              if(node==nod) then
                bis=.true.
                ee=0.5d0*(e0+e00)
                de=ee-e
              else if(node>nod) then
                de=-max(0.1d0,abs(e)/100.d0)
                e=emi
                tw=.false.
              else if(node<nod) then
                de=max(0.1d0,abs(e)/100.d0)
                e=ema
                tw=.false.
              endif
            else if(abs(q01).lt.abs(q001)) then
              e00=e0
              q001=q01
            else
              e=e-de
              de=-de/2
            endif
          endif
        else
          q1=q-1.d0
          if(q1*q001.lt.0.d0) then
            q01=q1
            e0=e
            e=0.5d0*(e00+e0)
            de=e-e0
            e=e0
          else
            q001=q1
            e00=e
            e=0.5d0*(e00+e0)
            de=e-e00
            e=e00 
          endif
        endif
      endif 
      if(e+de.ge.e2) then
        e=0.5d0*(e+e2)
      else if(e+de.le.e1) then
        e=0.5d0*(e+e1)
      else
        e=e+de
      endif
      fff=dabs(de/e)
      if(dabs(e).lt.1.d0) fff=dabs(de)
      if(fff.gt.tol.and.nit.lt.nitmax) goto 1
c     -----  normalize g -------
      fac=1.d0/dsqrt(q)
      do k=1,nre
        pnew(k)=pnew(k)*fac
        qnew(k)=qnew(k)*fac
      enddo
      do k=nre+1,nrad
        pnew(k)=0.d0
        qnew(k)=0.d0
      enddo
      if(key1==0) q=1.d0
      fi00=pnew(nrad)/r(nrad)
      h12=1.d0/12.d0/h
      if(irel0.eq.0) then
        dp=(l+1)/r(nrad)*pnew(nrad)+qnew(nrad)
        dfi00=dp/r(nrad)-pnew(nrad)/r(nrad)**2
        if(fi00/=0.d0) dl=dfi00/fi00*r(nrad)
        if(fi00==0.d0) dl=0.d0
      else if(irel0.eq.1) then
        work(nrad)=v(nrad)-(z+z)/r(nrad)-e
        alf=1.d0-work(nrad)/c2
        fac=dfloat(l)*(l+1.d0)/alf/r(nrad)**2
        dp=pnew(nrad)/r(nrad)+alf*qnew(nrad)
        if(key1.eq.1) dp=dp-f_exq(nrad)/c2
        dfi00=dp/r(nrad)-pnew(nrad)/r(nrad)**2
        if(fi00/=0.d0) dl=dfi00/fi00*r(nrad)
        if(fi00==0.d0) dl=0.d0
      else if(irel0.eq.2) then
        work(nrad)=v(nrad)-(z+z)/r(nrad)-e
        alf=1.d0-work(nrad)/c2
        dp=-akap*pnew(nrad)/r(nrad)+alf*qnew(nrad)
        if(key1.eq.1) dp=dp-f_exq(nrad)/c2
        dfi_big=dp/r(nrad)-pnew(nrad)/r(nrad)**2
        if(fi00/=0.d0) dl=dfi_big/fi00*r(nrad)
        if(fi00==0.d0) dl=0.d0
        dfi00=qnew(nrad)/r(nrad)
      endif
c     ------ Check of numerical accuracy ----------------
c     ------ by comparing numerical derivatives with the right parts ----
      acc=0.d0
      do ir=1,nrad
        if(ir==0) then
          p0=h12*(-25.d0*pnew(0)+48.d0*pnew(1)-36.d0*pnew(2)
     &      +16.d0*pnew(3)-3.d0*pnew(4))/dr(ir)
          q0=h12*(-25.d0*qnew(0)+48.d0*qnew(1)-36.d0*qnew(2)
     &      +16.d0*qnew(3)-3.d0*qnew(4))/dr(ir)
        else if(ir==1) then
          p0=h12*(-3.d0*pnew(0)-10.d0*pnew(1)+18.d0*pnew(2)
     &      -6.d0*pnew(3)+pnew(4))/dr(ir)
          q0=h12*(-3.d0*qnew(0)-108.d0*qnew(1)+18.d0*qnew(2)
     &      -6.d0*qnew(3)+qnew(4))/dr(ir)
        else if(ir==nrad-1) then
          p0=h12*(-pnew(nrad-4)+6.d0*pnew(nrad-3)-18.d0*pnew(nrad-2)
     &      +10.d0*pnew(nrad-1)+3.d0*pnew(nrad))/dr(ir)
          q0=h12*(-qnew(nrad-4)+6.d0*qnew(nrad-3)-18.d0*qnew(nrad-2)
     &      +10.d0*qnew(nrad-1)+3.d0*qnew(nrad))/dr(ir)
        else if(ir==nrad) then
          p0=h12*(3.d0*pnew(nrad-4)-16.d0*pnew(nrad-3)
     &      +36.d0*pnew(nrad-2)-48.d0*pnew(nrad-1)
     &      +25.d0*pnew(nrad))/dr(ir)
          q0=h12*(3.d0*qnew(nrad-4)-16.d0*qnew(nrad-3)
     &      +36.d0*qnew(nrad-2)-48.d0*qnew(nrad-1)
     &      +25.d0*qnew(nrad))/dr(ir)
        else
          p0=h12*(pnew(ir-2)+8.d0*(pnew(ir+1)-pnew(ir-1))-pnew(ir+2))
     &      /dr(ir)
          q0=h12*(qnew(ir-2)+8.d0*(qnew(ir+1)-qnew(ir-1))-qnew(ir+2))
     &      /dr(ir)
        endif
        if(irel0.eq.0) then
          p1=(l+1)/r(ir)*pnew(ir)+qnew(ir)
          q1=(v(ir)-(z+z)/r(ir)-e)*pnew(ir)-(l+1)/r(ir)*qnew(ir)
          if(key1.eq.1) q1=q1+f_ex(ir)
        else if(irel0.eq.1) then
          work(ir)=v(ir)-(z+z)/r(ir)-e
          alf=1.d0-work(ir)/c2
          p1=pnew(ir)/r(ir)+alf*qnew(ir)
          if(key1.eq.1) p1=p1-f_exq(ir)/c2
          fac=dfloat(l)*(l+1.d0)/alf/r(ir)**2
          q1=-qnew(ir)/r(ir)+(work(ir)+fac)*pnew(ir)
          if(key1.eq.1) q1=q1+f_ex(ir)
        else if(irel0.eq.2) then
          work(ir)=v(ir)-(z+z)/r(ir)-e
          alf=1.d0-work(ir)/c2
          p1=-akap*pnew(ir)/r(ir)+alf*qnew(ir)
          if(key1.eq.1) p1=p1-f_exq(ir)/c2
          q1=akap*qnew(ir)/r(ir)+work(ir)*pnew(ir)
          if(key1.eq.1) q1=q1+f_ex(ir)
        endif
        acc=max(acc,abs(p1-p0))
        acc=max(acc,abs(q1-q0))
      enddo
      end

      subroutine fit_zero(f,r)
      implicit none
      double precision, intent(inout) :: f(0:3),r(0:3)
      integer :: i
      double precision :: a(3),b(3),c(3),det3,det_all
      a=1.d0
      b=r(1:)
      do i=1,3
        c(i)=r(i)**2
      enddo
      det_all=det3(a,b,c)
      f(0)=det3(f(1),b,c)/det_all
      end



      subroutine get_nre(e,nre,l,rofi,v,z,nr)
      implicit none
      integer :: nre,nr,ir,l
      double precision :: v(0:nr),rofi(0:nr),fllp1,alf,bet,e,z,ex
      fllp1=l*(l+1)
      do ir=nr,3,-1
        alf=fllp1/rofi(ir)**2+v(ir)-2.d0*z/rofi(ir)-e
        if(alf.le.0.d0) then
          nre=ir
          goto 1
        endif
        bet=sqrt(alf)
c     if(bet*rofi(ir)<13.d0*log(10.d0)) then
c     nre=ir
c     goto 1
c     endif
        ex=exp(-bet*rofi(ir))
        if(ex.gt.1.d-10) then
          nre=ir
          goto 1
        endif
      enddo
 1    continue
      nre=(nre/2)*2
      end



      subroutine pn_qn(irel0,e,z,l,i,val,slo,v,r,pn,qn,f_exq,key1,
     $  akap,clight)

      implicit none
      integer, intent(in) :: irel0,l,i,key1
      double precision, intent(in) :: v,e,z,val,slo,r,f_exq,clight
      double precision, intent(out) :: pn,qn,akap
      integer :: kap
      double precision :: c2,work,alf
      c2=clight**2
      if(irel0.eq.0) then
        akap=dfloat(l+1)
        pn=val                  !! =P in V-19-3
        qn=slo-akap*val/r       !! =Q in V-19-3
      else if(irel0.eq.1) then
        akap=dfloat(l*(l+1))
        work=v-(z+z)/r-e
        alf=1.d0-work/c2
        pn=val                  !! =P in V-19-3
        if(key1==0) qn=(slo-val/r)/alf 
        if(key1==1) qn=(slo-val/r+f_exq/c2)/alf !in VI-87-28
      else if(irel0.eq.2) then
        kap=0
        if(i.lt.0) kap=l
        if(i.gt.0) kap=-l-1
        akap=kap
        work=v-(z+z)/r-e
        alf=1.d0-work/c2
        pn=val                  !! =P in V-19-3
        if(key1==0) qn=(slo+akap*val/r)/alf
        if(key1==1) qn=(slo+akap*val/r+f_exq/c2)/alf !87-17
      endif
      end



      SUBROUTINE RADSCH_b(e,z,a,b,l,m,r,dr,p,q,x_q,x_p,v,key1,akap,
     &  p1,q1,p0,q0,key2,nrad,nod2,q2,work,irel0,clight)
***************************************************************
*     The solution of the radial Schroedinger equation           *
*     *
*     / P'(r)=Q(r)+(l+1)*P(r)/r                                 *
*     \ Q'(r)=(V(r)-E)*P(r)-(l+1)*Q(r)/r  + Xp                  *
*     FI(r)=P(r)/r and FID(r)=FI'(r)                          *
*     *
*     or scalar-relativistic Dirac equation                      *
*     *
*     / P'(r)=A(r)*Q(r)+P(r)/r-Xq/c^2                           *
*     \ Q'(r)=(V(r)-E+C(r))*P(r)-Q(r)/r + Xp                    *
*     where                                                      *
*     A(r)=1+(E-V(r))/c^2                                     *
*     C(r)=l*(l+1)/r^2/A(r)                                   *
*     FI(r)=P(r)/r and FID(r)=FI'(r)                          *
*     and c is the light velocity                                *
***************************************************************
c     key1 = 1  Exchange is included
c     key2 = 1  Energy derivative finding
      IMPLICIT none
      integer, intent(in) :: irel0,nrad
      integer :: l,m,key1,key2,nod2,l1,i0,i,kc,ir,ii,k
      double precision :: v(0:nrad),P(0:nrad),Q(0:nrad),
     $  R(0:nrad),DR(0:nrad),
     &  x_p(0:nrad),x_q(0:nrad),q0(0:nrad),p0(0:nrad),
     &  work(0:nrad),pp(5),qq(5),e,z,a,b,akap,c2,z2,dfk2,r3,a3,
     &  b3,v3,xp3,yp3,r1,r2,a1,a2,b1,b2,v1,v2,xp1,xp2,yp1,yp2,
     &  ak1,am1,ak2,am2,ak3,am3,ak4,am4,ac2,ev3,alf,xq3,yq3,
     &  adr1,adr2,adr3,ev1,ev2,xq1,xq2,yq1,yq2,q2,p1,q1,p02,q02,
     &  adr,ev,fll,p_p,q_p,ppp,qqq,p_c,q_c,dqdall,clight
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
C     ****** NON-RELATIVISTIC VERSION (LOGARITHMIC SCALE) ******
      IF(irel0.EQ.0) THEN
C     ASYMPTOTIC BEHAVIOR ::
        I0=nrad
        p(i0)=p1
        q(i0)=q1
C     4-TH ORDER RUNGE-KUTTA LOOP ::
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
     &      (e+z2/r2-(9.D0*(v(I)+v(I+1))-v(I-1)-v(I+2))*6.25D-2)
          IF(I.EQ.nrad-1) V2=A2*(E+z2/r2-(v(I)+v(I+1))*0.5D0)
          V3=A3*(E+z2/r3-v(I))
          if(key1.eq.1) then
            xp1=xp3
            IF(I.NE.nrad-1) xp2=A2*
     &        (9.D0*(x_p(I)+x_p(I+1))-x_p(I-1)-x_p(I+2))*6.25D-2
            IF(I.EQ.nrad-1) xp2=A2*(x_p(I)+x_p(I+1))*0.5D0
            xp3=x_p(i)*a3
          endif
          if(key2.eq.1) then
            yp1=yp3
            IF(I.NE.nrad-1) yp2=A2*
     &        (9.D0*(p0(I)+p0(I+1))-p0(I-1)-p0(I+2))*6.25D-2
            IF(I.EQ.nrad-1) yp2=A2*(p0(I)+p0(I+1))*0.5D0
            yp3=p0(i)*a3
          endif
C     **************
          aK1= A1*Q(I+1)+B1*P(I+1)
          aM1=-V1*P(I+1)-B1*Q(I+1)
          if(key1.eq.1) am1=am1+xp1
          if(key2.eq.1) am1=am1-yp1
C     **************
          aK2= A2*(Q(I+1)-aM1*0.5D0)+B2*(P(I+1)-aK1*0.5D0)
          aM2=-V2*(P(I+1)-aK1*0.5D0)-B2*(Q(I+1)-aM1*0.5D0)
          if(key1.eq.1) am2=am2+xp2
          if(key2.eq.1) am2=am2-yp2
C     **************
          aK3= A2*(Q(I+1)-aM2*0.5D0)+B2*(P(I+1)-aK2*0.5D0)
          aM3=-V2*(P(I+1)-aK2*0.5D0)-B2*(Q(I+1)-aM2*0.5D0)
          if(key1.eq.1) am3=am3+xp2
          if(key2.eq.1) am3=am3-yp2
C     **************
          aK4= A3*(Q(I+1)-aM3)+B3*(P(I+1)-aK3)
          aM4=-V3*(P(I+1)-aK3)-B3*(Q(I+1)-aM3)
          if(key1.eq.1) am4=am4+xp3
          if(key2.eq.1) am4=am4-yp3
C     **************
          P(I)=P(I+1)-(aK1+aK4+2.D0*(aK2+aK3))/6.D0
          Q(I)=Q(I+1)-(aM1+aM4+2.D0*(aM2+aM3))/6.D0
        enddo
C     ****** RELATIVISTIC VERSION ******
      ELSEIF(irel0.EQ.1)THEN
C     ASyMPTOTIC BEHAVIOR ::
        ac2=akap/c2
        I0=nrad
        p(i0)=p1
        q(i0)=q1
C     4-TH ORDER RUNGE-KUTTA LOOP ::
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
          R2=B*(EXP(a*(I+0.5D0))-1.D0) !!! R(I+1/2)
          R3=R(I)
          adr1=adr3
          adr2=(r2+b)*a
          adr3=dr(i)*a
          B1=B3
          B2=adr2/R2
          B3=adr3/R3
          EV1=EV3
          IF(I.NE.nrad-1) eV2=e+z2/r2
     &      -(9.D0*(v(I)+v(I+1))-v(I-1)-v(I+2))*6.25D-2
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
     &        (9.D0*(x_p(I)+x_p(I+1))-x_p(I-1)-x_p(I+2))*6.25D-2
            IF(I.EQ.nrad-1) xp2=adr2*(x_p(I)+x_p(I+1))*0.5D0
            xp3=x_p(i)*adr3
            xq1=xq3
            IF(I.NE.nrad-1) xq2=adr2/c2*
     &        (9.D0*(x_q(I)+x_q(I+1))-x_q(I-1)-x_q(I+2))*6.25D-2
            IF(I.EQ.nrad-1) xq2=adr2/c2*(x_q(I)+x_q(I+1))*0.5D0
            xq3=x_q(i)*adr3/c2
          endif
          if(key2.eq.1) then
            yp1=yp3
            IF(I.NE.nrad-1) p02=
     &        (9.D0*(p0(I)+p0(I+1))-p0(I-1)-p0(I+2))*6.25D-2
            IF(I.EQ.nrad-1) p02=(p0(I)+p0(I+1))*0.5D0
            yp2=(1.d0+ac2/(alf*r2)**2)*p02*adr2
            yp3=(1.d0+ac2/(alf*r3)**2)*p0(i)*adr3
            yq1=yq3
            IF(I.NE.nrad-1) q02=
     &        (9.D0*(q0(I)+q0(I+1))-q0(I-1)-q0(I+2))*6.25D-2
            IF(I.EQ.nrad-1) q02=(q0(I)+q0(I+1))*0.5D0
            yq2=q02*adr2/c2
            yq3=q0(i)*adr3/c2
          endif
C     **************
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
C     **************
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
C     **************
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
C     **************
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
C     **************
          P(I)=P(I+1)-(aK1+aK4+2.D0*(aK2+aK3))/6.D0
          Q(I)=Q(I+1)-(aM1+aM4+2.D0*(aM2+aM3))/6.D0
        enddo
      ELSEIF(irel0.EQ.2)THEN
C     ASyMPTOTIC BEHAVIOR ::
        I0=nrad
        p(i0)=p1
        q(i0)=q1
C     4-TH ORDER RUNGE-KUTTA LOOP ::
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
          R2=B*(EXP(a*(I+0.5D0))-1.D0) !!! R(I+1/2)
          R3=R(I)
          adr1=adr3
          adr2=(r2+b)*a
          adr3=dr(i)*a
          B1=B3
          B2=-akap*adr2/R2
          B3=-akap*adr3/R3
          EV1=EV3
          IF(I.NE.nrad-1) eV2=e+z2/r2
     &      -(9.D0*(v(I)+v(I+1))-v(I-1)-v(I+2))*6.25D-2
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
     &        (9.D0*(x_p(I)+x_p(I+1))-x_p(I-1)-x_p(I+2))*6.25D-2
            IF(I.EQ.nrad-1) xp2=adr2*(x_p(I)+x_p(I+1))*0.5D0
            xp3=x_p(i)*adr3
            xq1=xq3
            IF(I.NE.nrad-1) xq2=adr2/c2*
     &        (9.D0*(x_q(I)+x_q(I+1))-x_q(I-1)-x_q(I+2))*6.25D-2
            IF(I.EQ.nrad-1) xq2=adr2/c2*(x_q(I)+x_q(I+1))*0.5D0
            xq3=x_q(i)*adr3/c2
          endif
          if(key2.eq.1) then
            yp1=yp3
            IF(I.NE.nrad-1) p02=
     &        (9.D0*(p0(I)+p0(I+1))-p0(I-1)-p0(I+2))*6.25D-2
            IF(I.EQ.nrad-1) p02=(p0(I)+p0(I+1))*0.5D0
            yp2=p02*adr2
            yp3=p0(i)*adr3
            yq1=yq3
            IF(I.NE.nrad-1) q02=
     &        (9.D0*(q0(I)+q0(I+1))-q0(I-1)-q0(I+2))*6.25D-2
            IF(I.EQ.nrad-1) q02=(q0(I)+q0(I+1))*0.5D0
            yq2=q02*adr2/c2
            yq3=q0(i)*adr3/c2
          endif
C     **************
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
C     **************
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
C     **************
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
C     **************
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
C     **************
          P(I)=P(I+1)-(aK1+aK4+2.D0*(aK2+aK3))/6.D0
          Q(I)=Q(I+1)-(aM1+aM4+2.D0*(aM2+aM3))/6.D0
        enddo
      ENDIF


      do i=1,5
        ii=nrad-i+1
        adr=a*dr(ii)
c     ------ See V-19 for details -------------------------------------------
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
c     if(key2.eq.2) qq(i)=qq(i)+2.d0*fll/alf**2/c2**2*p1(ii)
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
     &    -1274.d0*pp(2)+251.d0*pp(1))/720.d0
        q_p=q(i)-(1901.d0*qq(5)-2774.d0*qq(4)+2616.d0*qq(3)
     &    -1274.d0*qq(2)+251.d0*qq(1))/720.d0
        adr=a*dr(i-1)
c     ------ See V-19 for details -------------------------------------------
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
c     if(key2.eq.2) qqq=qqq+2.d0*fll/alf**2/c2**2*p1(i-1)
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
     &    +106.d0*pp(3)-19.d0*pp(2))/720.d0
        p(i-1)=(475.d0*p_c+27.d0*p_p)/502.d0
        q_c=q(i)-(251.d0*qqq+646.d0*qq(5)-264.d0*qq(4)
     &    +106.d0*qq(3)-19.d0*qq(2))/720.d0
        q(i-1)=(475.d0*q_c+27.d0*q_p)/502.d0
        do k=1,4
          pp(k)=pp(k+1)
          qq(k)=qq(k+1)
        enddo
        adr=a*dr(i-1)
c     ------ See V-19 for details -------------------------------------------
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
c     if(key2.eq.2) qq(5)=qq(5)+2.d0*fll/alf**2/c2**2*p1(i-1)
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
 2    continue
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
      


      SUBROUTINE RADSCH(e,z,a,b,l,m,r,dr,p,q,x_q,x_p,v,key1,akap,key,
     &  p11,q11,p0,q0,key2,nrad,nod1,qq1,work,p1,irel0,clight)
***************************************************************
*     The solution of the radial Schroedinger equation           *
*     *
*     / P'(r)=Q(r)+(l+1)*P(r)/r                                 *
*     \ Q'(r)=(V(r)-E)*P(r)-(l+1)*Q(r)/r  + Xp                  *
*     FI(r)=P(r)/r and FID(r)=FI'(r)                          *
*     *
*     or scalar-relativistic Dirac equation                      *
*     *
*     / P'(r)=A(r)*Q(r)+P(r)/r-Xq/c^2                           *
*     \ Q'(r)=(V(r)-E+C(r))*P(r)-Q(r)/r + Xp                    *
*     where                                                      *
*     A(r)=1+(E-V(r))/c^2                                     *
*     C(r)=l*(l+1)/r^2/A(r)                                   *
*     FI(r)=P(r)/r and FID(r)=FI'(r)                          *
*     and c is the light velocity                                *
***************************************************************
c     key  = 1  Old values P and Q are used to start the integration
c     key1 = 1  Exchange is included
c     key2 = 1  Energy derivative finding
c     key2 = 2  Second energy derivative finding
      IMPLICIT none
      integer, intent(in) :: nrad,key,key1,key2,irel0
      double precision, intent(in) :: p1(0:nrad)
      integer :: l,m,nod1,l1,m_rk,mn,i,i0,i00,k,ir
      double precision :: v(0:nrad),P(0:nrad),Q(0:nrad),
     $  R(0:nrad),DR(0:nrad),
     &  x_p(0:nrad),x_q(0:nrad),q0(0:nrad),p0(0:nrad),
     &  work(0:nrad),arg(10),pp(5),qq(5),e,z,a,b,akap,c2,z2,
     &  dfk2,r3,a3,
     &  b3,v3,xp3,yp3,r1,r2,a1,a2,b1,b2,v1,v2,xp1,xp2,yp1,yp2,
     &  ak1,am1,ak2,am2,ak3,am3,ak4,am4,ac2,ev3,alf,xq3,yq3,
     &  adr1,adr2,adr3,ev1,ev2,xq1,xq2,yq1,yq2,p11,q11,t0,vvv,
     &  ppp,ql,p02,q02,z2c2,alfa,a0,b0,adr,ev,fll,p_p,q_p,qqq,
     &  p_c,q_c,qq1,dqdall,clight
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
      m_rk=5                    !! number of Runge-Kutta steps
      mn=6
      do i=1,mn
        arg(i)=dfloat(i-1)
      enddo
C     ****** NON-RELATIVISTIC VERSION (LOGARITHMIC SCALE) ******
      IF(irel0.EQ.0) THEN
C     ASSYMPTOTIC BEHAVIOR ::
        P(0)=0.D0
        Q(0)=0.D0
        I0=1
        if(key.eq.0) then
 1        P(I0)=R(I0)**L1*(1.D0-z2/(2.D0*L1)*R(I0))
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
c     IF(Z.LE.1.D-05)Q(I0)=-(vnuc+EN)*R(I0)**(L1+1)/(2.D0*L1+1)
        else
          p(i0)=p11
          q(i0)=q11
        endif
C     4-TH ORDER RUNGE-KUTTA LOOP ::
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
     &      (e+z2/r2-(9.D0*(v(I)+v(I+1))-v(I-1)-v(I+2))*6.25D-2)
          IF(I.EQ.nrad-1) V2=A2*(E+z2/r2-(v(I)+v(I+1))*0.5D0)
          v2=a2*(e+z2/r2-vvv)
          V3=A3*(E+z2/r3-v(I+1))
          if(key1.eq.1) then
            xp1=xp3
            IF(I.NE.nrad-1) xp2=A2*
     &        (9.D0*(x_p(I)+x_p(I+1))-x_p(I-1)-x_p(I+2))*6.25D-2
            IF(I.EQ.nrad-1) xp2=A2*(x_p(I)+x_p(I+1))*0.5D0
            xp3=x_p(i+1)*a3
          endif
          if(key2.ge.1) then
            yp1=yp3
	    call fint(p0(i00),arg,mn,t0,ppp)
            IF(I.NE.nrad-1) yp2=A2*dfk2*
     &        (9.D0*(p0(I)+p0(I+1))-p0(I-1)-p0(I+2))*6.25D-2
            IF(I.EQ.nrad-1) yp2=dfk2*A2*(p0(I)+p0(I+1))*0.5D0
            yp2=a2*dfk2*ppp
            yp3=dfk2*p0(i+1)*a3
          endif
C     **************
          aK1= A1*Q(I)+B1*P(I)
          aM1=-V1*P(I)-B1*Q(I)
          if(key1.eq.1) am1=am1+xp1
          if(key2.ge.1) am1=am1-yp1
C     **************
          aK2= A2*(Q(I)+aM1*0.5D0)+B2*(P(I)+aK1*0.5D0)
          aM2=-V2*(P(I)+aK1*0.5D0)-B2*(Q(I)+aM1*0.5D0)
          if(key1.eq.1) am2=am2+xp2
          if(key2.ge.1) am2=am2-yp2
C     **************
          aK3= A2*(Q(I)+aM2*0.5D0)+B2*(P(I)+aK2*0.5D0)
          aM3=-V2*(P(I)+aK2*0.5D0)-B2*(Q(I)+aM2*0.5D0)
          if(key1.eq.1) am3=am3+xp2
          if(key2.ge.1) am3=am3-yp2
C     **************
          aK4= A3*(Q(I)+aM3)+B3*(P(I)+aK3)
          aM4=-V3*(P(I)+aK3)-B3*(Q(I)+aM3)
          if(key1.eq.1) am4=am4+xp3
          if(key2.ge.1) am4=am4-yp3
C     **************
          P(I+1)=P(I)+(aK1+aK4+2.D0*(aK2+aK3))/6.D0
          Q(I+1)=Q(I)+(aM1+aM4+2.D0*(aM2+aM3))/6.D0
        enddo
C     ****** RELATIVISTIC VERSION ******
      ELSEIF(irel0.EQ.1)THEN
C     ASSIMPTOTIC BEHAVIOR ::
        ac2=akap/c2
        QL=SQRT(akap+1.d0-z2**2/C2)
        P(0)=0.D0
        Q(0)=0.D0
        I0=1
        if(key.eq.0) then
 21       P(I0)=R(I0)**QL
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
C     4-TH ORDER RUNGE-KUTTA LOOP ::
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
          R2=B*(EXP(a*(I+0.5D0))-1.D0) !!! R(I+1/2)
          R3=R(I+1)
          adr1=adr3
          adr2=(r2+b)*a
          adr3=dr(i+1)*a
          B1=B3
          B2=adr2/R2
          B3=adr3/R3
          EV1=EV3
          IF(I.NE.nrad-1) eV2=e+z2/r2
     &      -(9.D0*(v(I)+v(I+1))-v(I-1)-v(I+2))*6.25D-2
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
     &        (9.D0*(x_p(I)+x_p(I+1))-x_p(I-1)-x_p(I+2))*6.25D-2
            IF(I.EQ.nrad-1) xp2=adr2*(x_p(I)+x_p(I+1))*0.5D0
            xp3=x_p(i+1)*adr3
            xq1=xq3
            IF(I.NE.nrad-1) xq2=adr2/c2*
     &        (9.D0*(x_q(I)+x_q(I+1))-x_q(I-1)-x_q(I+2))*6.25D-2
            IF(I.EQ.nrad-1) xq2=adr2/c2*(x_q(I)+x_q(I+1))*0.5D0
            xq3=x_q(i+1)*adr3/c2
          endif
          if(key2.ge.1) then
            yp1=yp3
            IF(I.NE.nrad-1) p02=dfk2*
     &        (9.D0*(p0(I)+p0(I+1))-p0(I-1)-p0(I+2))*6.25D-2
            IF(I.EQ.nrad-1) p02=dfk2*(p0(I)+p0(I+1))*0.5D0
            yp2=(1.d0+ac2/(alf*r2)**2)*p02*adr2
            yp3=dfk2*(1.d0+ac2/(alf*r3)**2)*p0(i+1)*adr3
            yq1=yq3
            IF(I.NE.nrad-1) q02=
     &        (9.D0*(q0(I)+q0(I+1))-q0(I-1)-q0(I+2))*6.25D-2
            IF(I.EQ.nrad-1) q02=(q0(I)+q0(I+1))*0.5D0
            yq2=dfk2*q02*adr2/c2
            yq3=dfk2*q0(i+1)*adr3/c2
          endif
C     **************
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
C     **************
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
C     **************
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
C     **************
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
C     **************
          P(I+1)=P(I)+(aK1+aK4+2.D0*(aK2+aK3))/6.D0
          Q(I+1)=Q(I)+(aM1+aM4+2.D0*(aM2+aM3))/6.D0
        enddo
      ELSEIF(irel0.EQ.2)THEN
C     ASSIMPTOTIC BEHAVIOR ::
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
 31       P(I0)=a0*R(I0)**alfa
c     IF(P(I0).EQ.0.D0)THEN
c     Q(I0)=0.D0
c     I0=I0+1
c     GOTO 31
c     ENDIF
          Q(I0)=b0*r(i0)**alfa
        else
          p(i0)=p11
          q(i0)=q11
        endif
C     4-TH ORDER RUNGE-KUTTA LOOP ::
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
          R2=B*(EXP(a*(I+0.5D0))-1.D0) !!! R(I+1/2)
          R3=R(I+1)
          adr1=adr3
          adr2=(r2+b)*a
          adr3=dr(i+1)*a
          B1=B3
          B2=-akap*adr2/R2
          B3=-akap*adr3/R3
          EV1=EV3
          IF(I.NE.nrad-1) eV2=e+z2/r2
     &      -(9.D0*(v(I)+v(I+1))-v(I-1)-v(I+2))*6.25D-2
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
     &        (9.D0*(x_p(I)+x_p(I+1))-x_p(I-1)-x_p(I+2))*6.25D-2
            IF(I.EQ.nrad-1) xp2=adr2*(x_p(I)+x_p(I+1))*0.5D0
            xp3=x_p(i+1)*adr3
            xq1=xq3
            IF(I.NE.nrad-1) xq2=adr2/c2*
     &        (9.D0*(x_q(I)+x_q(I+1))-x_q(I-1)-x_q(I+2))*6.25D-2
            IF(I.EQ.nrad-1) xq2=adr2/c2*(x_q(I)+x_q(I+1))*0.5D0
            xq3=x_q(i+1)*adr3/c2
          endif
          if(key2.ge.1) then
            yp1=yp3
            IF(I.NE.nrad-1) p02=
     &        (9.D0*(p0(I)+p0(I+1))-p0(I-1)-p0(I+2))*6.25D-2
            IF(I.EQ.nrad-1) p02=(p0(I)+p0(I+1))*0.5D0
            yp2=dfk2*p02*adr2
            yp3=dfk2*p0(i+1)*adr3
            yq1=yq3
            IF(I.NE.nrad-1) q02=
     &        (9.D0*(q0(I)+q0(I+1))-q0(I-1)-q0(I+2))*6.25D-2
            IF(I.EQ.nrad-1) q02=(q0(I)+q0(I+1))*0.5D0
            yq2=dfk2*q02*adr2/c2
            yq3=dfk2*q0(i+1)*adr3/c2
          endif
C     **************
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
C     **************
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
C     **************
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
C     **************
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
C     **************
          P(I+1)=P(I)+(aK1+aK4+2.D0*(aK2+aK3))/6.D0
          Q(I+1)=Q(I)+(aM1+aM4+2.D0*(aM2+aM3))/6.D0
        enddo
      ENDIF

      
      do ir=1,5
        i=m_rk+ir-5
        adr=a*dr(i)
c     ------ See V-19 for details -------------------------------------------
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
     &    -1274.d0*pp(2)+251.d0*pp(1))/720.d0
        q_p=q(i)+(1901.d0*qq(5)-2774.d0*qq(4)+2616.d0*qq(3)
     &    -1274.d0*qq(2)+251.d0*qq(1))/720.d0
        adr=a*dr(i+1)
c     ------ See V-19 for details -------------------------------------------
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
     &    +106.d0*pp(3)-19.d0*pp(2))/720.d0
        p(i+1)=(475.d0*p_c+27.d0*p_p)/502.d0
        q_c=q(i)+(251.d0*qqq+646.d0*qq(5)-264.d0*qq(4)
     &    +106.d0*qq(3)-19.d0*qq(2))/720.d0
        q(i+1)=(475.d0*q_c+27.d0*q_p)/502.d0
        do k=1,4
          pp(k)=pp(k+1)
          qq(k)=qq(k+1)
        enddo
        adr=a*dr(i+1)
c     ------ See V-19 for details -------------------------------------------
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

      subroutine inhomo(ph_in,pi_in,qh_in,qi_in,ph,qh,ppi,qqi,node,
     &  q,nr,kc,dr,h,irel,clight)
      implicit none
      integer, intent(in) :: nr,kc,irel
      integer, intent(out) :: node
      double precision, intent(in) :: ph_in,pi_in,qh_in,qi_in,
     $  ppi(0:nr),qqi(0:nr),
     &  dr(0:nr),h,clight
      double precision, intent(inout) :: ph(0:nr),qh(0:nr)
      double precision, intent(out) :: q
      integer :: ir
      double precision :: ph_out,qh_out,pi_out,qi_out,dt,pp,qq,aa,bb,c2
      double precision, allocatable :: work(:)
      double precision, external :: dqdall      
      ph_out=ph(kc)
      qh_out=qh(kc)
      pi_out=ppi(kc)
      qi_out=qqi(kc)
      dt=qh_out*ph_in-ph_out*qh_in
      pp=pi_in-pi_out
      qq=qi_in-qi_out
      aa=(ph_in*qq-pp*qh_in)/dt
      bb=(ph_out*qq-pp*qh_out)/dt
      do ir=0,kc
        ph(ir)=ppi(ir)+aa*ph(ir)
        qh(ir)=qqi(ir)+aa*qh(ir)
      enddo
      do ir=kc+1,nr
        ph(ir)=ppi(ir)+bb*ph(ir)
        qh(ir)=qqi(ir)+bb*qh(ir)
      enddo
      allocate(work(0:nr))
      do ir=0,nr
        work(ir)=ph(ir)**2*dr(ir)
      enddo
      if(irel.ge.1) then
        c2=clight**2
        do ir=0,nr
          work(ir)=work(ir)+qh(ir)**2*dr(ir)/c2
        enddo
      endif
      q=dqdall(h,work(0),nr)
      deallocate(work)
      node=0
      do ir=1,nr
        if(ph(ir-1)*ph(ir).lt.0.d0) node=node+1
      enddo
      end      


      subroutine rad_hf_check(pn,qn,nre,node,nr)
      implicit none
      integer, intent(in) :: nr
      integer, intent(inout) :: nre,node
      double precision, intent(inout) :: pn(0:nr),qn(0:nr)
      integer :: i,ii
      double precision :: a
      integer, allocatable :: nn(:)
      double precision, allocatable :: p(:)
      allocate(p(node+1))
      allocate(nn(node+1))
      ii=1
      nn=0
      p=0.d0
      a=pn(1)
      do i=2,nre
        if(a*pn(i)<0.d0) then
          nn(ii)=i-1
          ii=ii+1
        endif
        a=pn(i)
        p(ii)=max(p(ii),abs(a))
      enddo
      nn(ii)=nre
      do i=2,ii
        if(p(i)/p(i-1)<1.d-1) then
          node=i-2
          nre=nn(i-1)
          goto 1
        endif
      enddo
 1    deallocate(p,nn)
      pn(nre+1:nr)=0.d0
      qn(nre+1:nr)=0.d0
      if(abs(pn(nre))>1.d-7.and.nre<nr) nre=nre+1
      end      



      subroutine fint(f,x,n,t,p)
c     interpolation of function f(x) in point t.
      implicit none
      integer, intent(in) :: n
      double precision, intent(in) :: f(n),x(n),t
      double precision, intent(out) :: p
      integer :: i,j
      double precision :: s
      p=0.d0
      do i=1,n
        s=f(i)
        do j=1,i-1
          s=s*(t-x(j))/(x(i)-x(j))
        enddo
        do j=i+1,n
          s=s*(t-x(j))/(x(i)-x(j))
        enddo
        p=p+s
      enddo
      end      



      subroutine print_radial_function(ispin)

      use comwann_mod
      implicit none
      integer, intent(in) :: ispin
      integer :: isort,ir,ib,iatom
      character*100 :: ch_all, ch_ie, ch_li, ch_atom, ch_spin

      if (maswrk) then      
        do iatom=1, natom
          isort=is(iatom)          
          if (iatom .lt. 10) then
            write(ch_atom, '(i1)') iatom
          elseif (iatom .lt. 100) then
            write(ch_atom, '(i2)') iatom
          elseif (iatom .lt. 1000) then
            write(ch_atom, '(i3)') iatom
          endif
          
          if (nspin .eq. 1) then          
            ch_all
     $        =adjustl(trim(ch_atom))
          else
            write(ch_spin, '(i1)') ispin
            ch_all
     $        =adjustl(trim(ch_atom))//'_'
     $        //adjustl(trim(ch_spin))
          endif
          
          open(unit=181, file='radial_function_'
     $      //adjustl(trim(ch_all))//'.dat')
          do ir=1, nrad(isort)        
            write(181, '(1000(e20.12,2x))') r(ir,isort),
     $        (bound_radfun(ir,1,ib,iatom,ispin), ib=1,
     $        nrel*lmb(isort)+1)
          enddo
          close(181)
        enddo
      endif
      end
      
      
      
      subroutine t_t1_x_all

      use comwann_mod
      integer :: isort, ispin
      double precision ::
     $  t1_x_h(0:maxb+maxb,maxlfun,maxlfun,nrel*maxb+1),
     $  t_x_h(nrel*maxb+1,nrel*maxl_cor+1,0:maxb+maxl_cor)
      
      t_x=0.0d0
      t1_x=0.0d0      
      
      do isort=1,nsort
        do ispin=1,nspin
          call t_t1_x(isort,ispin,1,t_x_h,t1_x_h)
          t1_x(:,:,:,:,isort,ispin)=t1_x_h
          t_x(:,:,:,isort,ispin)=t_x_h
        enddo
      enddo
      
      end


      subroutine t_t1_x(isort,ispin,key,t_x_h,t1_x_h)
c     key = 0 : for core levels
c     key = 1 : for valence levels

!     !!!!! complex*16 :: g_loc_0(:,:,:,:), 	allocate(g_loc_0(maxel,maxel,natom,nspin_0)) -> restart
!     !!!!! integer :: iat_1(:), iat_1(nsort)
!     !!!!! integer :: irel_core :: integer::  !!!!! 
!     !!!!! integer :: lc_max(:), lc_max(nsort)
!     !!!!! integer :: lm_isz(:,:) lm_isz(maxel, nsort)
!     !!!!  integer :: ncor(:) ncor(nsort)
!     !!!!! integer :: nrel_core:  !!!!!
!     !!!!! integer :: maxl_cor
!     !!!!! double precision :: t1_x(:,:,:,:)       allocate(t1_x(0:maxb+maxb,maxlfun,maxlfun,nrel*maxb+1))

!     !!!!! double precision :: t_x(:,:,:)  if(maxl_cor>=0) allocate(t_x(nrel*maxb+1,nrel*maxl_cor+1,0:maxb+maxl_cor))      
      use comwann_mod
      implicit none
      integer, intent(in) :: key,isort,ispin
      double precision, intent(out) ::
     $  t1_x_h(0:maxb+maxb,maxlfun,maxlfun,nrel*maxb+1),
     $  t_x_h(nrel*maxb+1,nrel*maxl_cor+1,0:maxb+maxl_cor)
      integer :: nlm_c,lm1,li1,lm,li,ii,iii,indx,
     &  icg1,icg2,icg3,icg4,icg,lm22,l22,km1,lf1,l1,i1,mj1,
     &  km2,lf2,li2,l2,i2,mj2,ii1,iii1,indx1,icgg,lm20,ic,l,
     &  lget,i,mj,lm2,iatom,lmc,lic,lc,mjc,nlm_all
      double precision :: anl,df22,fac,facj,cgg,pi8
      complex*16, allocatable :: gloc(:,:,:)
      allocate(gloc(maxel,maxel,natom))
      gloc=g_loc_0(:,:,:,ispin)
      pi8=8.d0*pi
      nlm_c=nrel_core*(lc_max(isort)+1)**2
      nlm_all=0
      t1_x_h=0.0d0
      t_x_h=0.0d0
      if(key==0) nlm_all=nlm_c
      if(key==1) nlm_all=nrel*(lmb(isort)+1)**2
      if(ncor(isort)/=0) then
        t_x_h=0.d0
        do lmc=1,nlm_c
          if(irel_core.le.1) then
            lic=lget(lmc)+1
          else if(irel_core.eq.2) then
            call getlimj(lmc,lc,ic,mjc,lic,0)
          endif
          do lm=1,nlm_all
            anl=0.0d0
            if(irel_core.le.1) then
              li=lget(lm)+1
              anl=dfloat(li+li-1)
            else if(irel_core.eq.2) then
              call getlimj(lm,l,i,mj,li,0)
              anl=dfloat(((li+1)/2)*2)
            endif
            if(irel_core.le.1) then
              ii = max0(lm,lmc)
              iii = min0(lm,lmc)
              indx = (ii*(ii-1))/2 + iii
              icg1 = indxcg(indx)
              icg2 = indxcg(indx+1) - 1 
            else
              icg1=indxcgr(lm,lmc,1,0)
              icg2=indxcgr(lm,lmc,2,0)
            endif
            do icg = icg1, icg2
              if(irel_core.le.1) then
                lm22=jcg(icg)
                cgg=cg(icg)**2
              else
                lm22=jcgr(icg)
                cgg=dconjg(cgr(icg))*cgr(icg)
              endif
              l22=lget(lm22)
              if(key==0.and.l22>2*lc_max(isort)) cycle
              if(key==1.and.l22>2*lmb(isort)) cycle
              df22=-dfloat(l22+l22+1)*anl/pi8
              t_x_h(li,lic,l22)=t_x_h(li,lic,l22)+cgg/df22
            enddo
          enddo                 !! over lm - core lm
        enddo                   !! over lmc
      endif
      t1_x_h=0.d0
      iatom=iat_1(isort)
      do km1=1,lfunm(isort) 
        lm1=lm_isz(km1,isort)
        lf1=lf_isz(km1,isort)
        if(irel.le.1) then
          li1=lget(lm1)+1
        else if(irel.eq.2) then
          call getlimj(lm1,l1,i1,mj1,li1,0)
        endif
        do lm=1,nlm_all
          anl=0.0d0
          if(irel_core.le.1) then
            li=lget(lm)+1
            anl=dfloat(li+li-1)
          else if(irel_core.eq.2) then
            call getlimj(lm,l,i,mj,li,0)
            anl=dfloat(((li+1)/2)*2)
          endif
          if(irel.le.1) then
            ii = max0(lm,lm1)
            iii = min0(lm,lm1)
            indx = (ii*(ii-1))/2 + iii
            icg1 = indxcg(indx)
            icg2 = indxcg(indx+1) - 1 
          else
            icg1=indxcgr(lm,lm1,1,0)
            icg2=indxcgr(lm,lm1,2,0)
          endif
          do icg = icg1, icg2
            if(irel.le.1) then
              lm22=jcg(icg)
              cgg=cg(icg)**2
            else
              lm22=jcgr(icg)
              cgg=dconjg(cgr(icg))*cgr(icg)
            endif
            l22=lget(lm22)
            df22=-dfloat(l22+l22+1)*anl/pi8
            do km2=1,lfunm(isort)
              lm2=lm_isz(km2,isort)
              lf2=lf_isz(km2,isort)
              if(irel.le.1) then
                li2=lget(lm2)+1
              else if(irel.eq.2) then
                call getlimj(lm2,l2,i2,mj2,li2,0)
              endif
              if(irel.le.1) then
                ii1 = max0(lm2,lm)
                iii1 = min0(lm2,lm)
                indx1 = (ii1*(ii1-1))/2 + iii1
                icg3 = indxcg(indx1)
                icg4 = indxcg(indx1+1) - 1 
              else
                icg3=indxcgr(lm2,lm,1,0)
                icg4=indxcgr(lm2,lm,2,0)
              endif
              do icgg = icg3, icg4
                if(irel.le.1) then
                  lm20 = jcg(icgg)
                else
                  lm20=jcgr(icgg)
                endif
                if(lm20.eq.lm22) then
                  if(irel.le.1) then
                    fac=cg(icg)*cg(icgg)/df22
                  else
                    fac=dreal(cgr(icg)*cgr(icgg)/df22)
                    facj=dimag(cgr(icg)*cgr(icgg)/df22)
                  endif
                  if(irel.le.1) then
                    t1_x_h(l22,lf2,lf1,li)=t1_x_h(l22,lf2,lf1,li)
     &                +dreal(gloc(km2,km1,iatom))*fac
                  else if(irel.eq.2) then
                    t1_x_h(l22,lf2,lf1,li)=t1_x_h(l22,lf2,lf1,li)
     &                +dreal(gloc(km2,km1,iatom)*dcmplx(fac,facj))
                  endif
                endif
              enddo
            enddo               !! over km2
          enddo
        enddo                   !! over lm - core lm
      enddo                     !! over km1
      deallocate(gloc)
      end



      subroutine f_ex_new(f_ex,f_exq,isort,pp,qq,ispin,li,nw)
!     !!!!! integer :: maxmtcor 
!     !!!!! integer :: ncormax
!     !!!! integer :: indcor(:,:,:) 	allocate(indcor(ncormax,nspin_0,nsort))
!     !!!! integer :: li_cor(:,:)  	allocate(li_cor(ncormax,nsort))
!     !!!!! integer :: max_rl 
!     !!!!! integer :: nspin_1
!     !!!!! double precision :: pcor(:) pcor(maxmtcor)-> restart
!     !!!!! double precision :: qcor(:) qcor(maxmtcor)-> restart
!     !!!!! double precision :: rl(:,:,:)       allocate(rl(0:maxnrad,0:max_rl,nsort))
      
      use comwann_mod
      implicit none
      integer, intent(in) :: isort,ispin,li,nw
      double precision, intent(in) :: pp(0:nw),qq(0:nw)
      double precision, intent(out) :: f_ex(0:nw),f_exq(0:nw)
      integer :: li1,lf1,l1,lf2,li2,l2,l,ic1,mt1,ir,mt2,ll2,isp
      double precision :: wh1
      double precision, allocatable :: pval(:),qval(:),pp1(:),f(:),g(:),
     $  qq1(:),tmp(:),work(:)
      isp=min(ispin,nspin_1)
      allocate(tmp(0:nw),pp1(0:nw),f(0:nw),g(0:nw),qval(0:nw),qq1(0:nw),
     &  work(0:nw),pval(0:nw))
      f_ex(0:nw)=0.d0
      f_exq(0:nw)=0.d0
      if(irel_core.le.1) then
        l=li-1
      else
        l=li/2
      endif
c     ------ core contribution --------------
      do ic1=1,ncor(isort)
        li1=li_cor(ic1,isort)
        if(irel_core.le.1) then
          l1=li1-1
        else
          l1=li1/2
        endif
        mt1=indcor(ic1,ispin,isort)
        call dcopy(nrad(isort)+1,pcor(mt1),1,pp1(0),1)
        if(irel.ge.1) call dcopy(nrad(isort)+1,qcor(mt1),1,qq1(0),1)
        do l2=0,min(l+l1,max_rl)
          call var_int(dr(0,isort),pp(0),qq(0),pp1(0),qq1(0),f(0),g(0),
     &      tmp(0),r(0,isort),rl(0,l2,isort),work(0),
     &      nrad(isort)+1,h(isort),clight,irel)
          wh1=t_x(li,li1,l2,isort,ispin)
          call dcopy(nrad(isort)+1,tmp(0),1,f(0),1)
          call dscal(nrad(isort)+1,wh1,f(0),1)
          do ir=0,nrad(isort)
            f_ex(ir)=f_ex(ir)+f(ir)*pp1(ir)
          enddo
          if(irel_core.ge.1) then
            do ir=0,nrad(isort)
              f_exq(ir)=f_exq(ir)+f(ir)*qq1(ir)
            enddo
          endif
        enddo                   !! over l2
      enddo                     !! over ic1
c     ------ Valence contribution ----------------------------------------
      do lf2=1,lfun(isort)
        li2=li_val(lf2,isort)
        if(irel.le.1) then
          ll2=li2-1
        else
          ll2=li2/2
        endif
        mt2=ind_wf(lf2,isort)
        do ir=0,nrad(isort)
          pval(ir)=gfun(mt2+ir,isp)*r(ir,isort)
        enddo
        if(irel.ge.1) then
          do ir=0,nrad(isort)
            qval(ir)=gfund(mt2+ir,isp)*r(ir,isort)
          enddo
        endif
        do l2=0,min(l+ll2,max_rl)
          call var_int(dr(0,isort),pval(0),qval(0),pp(0),qq(0),
     &      f(0),g(0),tmp(0),r(0,isort),rl(0,l2,isort),
     &      work(0),nrad(isort)+1,h(isort),clight,irel)
          do lf1=1,lfun(isort)
            mt1=ind_wf(lf1,isort)
            do ir=0,nrad(isort)
              f(ir)=gfun(mt1+ir,isp)*r(ir,isort)
              if(irel.ge.1) g(ir)=gfund(mt1+ir,isp)*r(ir,isort)
            enddo
            do ir=0,nrad(isort)
              f_ex(ir)=f_ex(ir)
     $          +tmp(ir)*t1_x(l2,lf2,lf1,li,isort,ispin)*f(ir)
            enddo
            if(irel>=1) then
              do ir=0,nrad(isort)
                f_exq(ir)=f_exq(ir)
     $            +tmp(ir)*t1_x(l2,lf2,lf1,li,isort,ispin)*g(ir)
              enddo
            endif
          enddo                 !! over lf1
        enddo                   !! over l2
      enddo                     !! over lf2
      deallocate(tmp,pp1,f,g,qval,qq1,work,pval)
      end
      


      subroutine var_int(ww,p,q,p0,q0,f,g,fg,r,rl,work,nr,h,clight,irel)
      implicit none
      integer :: nr,ir,irel
      double precision ::  ww(nr),p(nr),q(nr),p0(nr),q0(nr),
     $  f(nr),g(nr),fg(nr),rl(nr),r(nr),work(nr),c2,h,clight
      c2=clight**2
      do ir=1,nr
        work(ir)=p(ir)*p0(ir)
      enddo
      if(irel>=1) then
        do ir=1,nr
          work(ir)=work(ir)+q(ir)*q0(ir)/c2
        enddo
      endif
      do ir=1,nr
        work(ir)=ww(ir)*work(ir)*h
      enddo
      do ir=1,nr
        f(ir)=work(ir)*rl(ir)
      enddo
      call intalt(nr,f,g)
      fg(1)=0.d0
      do ir=2,nr
        fg(ir)=g(ir)/rl(ir)/r(ir)
      enddo
      f(1)=0.d0
      do ir=2,nr
        f(ir)=work(ir)/r(ir)/rl(ir)
      enddo
      call intalt(nr,f,g)
      do ir=1,nr
        fg(ir)=fg(ir)+(g(nr)-g(ir))*rl(ir)
      enddo
      end
      


      


      SUBROUTINE intalt(N,f,g)
      IMPLICIT none
      integer :: n,i
      double precision :: F(n),g(n),h12,h3,sum,f1,f2,f3
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
