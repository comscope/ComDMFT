      subroutine rad_dif(e,z,l,i,v,dv,r,nrad,p0,q0,p1,q1,p2,q2,dp0,
     &	               dq0,dp1,dq1,dp2,dq2,d2p0,d2q0,d2p1,d2q1,d2p2,
     &                   d2q2,kda)
c	The derivatives of radial solutions for given e
c	Input:
c	    p0,q0,p1,q1,p2,q2
c	Output:
c	    dp0,dq0,dp1,dq1,dp2,dq2
c	    d2p0,d2q0,d2p1,d2q1,d2p2,d2q2
	use manager_mod
	use units_mod
      implicit none
	integer :: l,i,nrad,kda,ir,kap
      real*8 :: v(0:nrad),dv(0:nrad),r(0:nrad),p0(0:nrad),q0(0:nrad),
     &          p1(0:nrad),q1(0:nrad),p2(0:nrad),q2(0:nrad),
     &          dp0(0:nrad),dq0(0:nrad),dp1(0:nrad),dq1(0:nrad),
     &          dp2(0:nrad),dq2(0:nrad),d2p0(0:nrad),d2q0(0:nrad),
     &          d2p1(0:nrad),d2q1(0:nrad),d2p2(0:nrad),d2q2(0:nrad),e,z,
     &          c2,z2,akap,r0,rr2,ve,dve,alf,dalf,a
      c2=clight**2
	z2=z+z
	if(irel.eq.0) then
	  akap=dfloat(l+1)
      else if(irel.eq.1) then
        akap=dfloat(l*(l+1))
      else if(irel.eq.2) then
        kap=0
        if(i.lt.0) kap=l
        if(i.gt.0) kap=-l-1
        akap=dfloat(kap)
	endif
	if(irel.eq.0) then
	  do ir=1,nrad
	    r0=r(ir)
	    rr2=r0*r0
	    ve=v(ir)-z2/r0-e
	    dve=dv(ir)+z2/rr2
	    dp0(ir)=q0(ir)+akap*p0(ir)/r0
	    dq0(ir)=ve*p0(ir)-akap*q0(ir)/r0
	    dp1(ir)=q1(ir)+akap*p1(ir)/r0
	    dq1(ir)=ve*p1(ir)-akap*q1(ir)/r0-p0(ir)
	    d2p0(ir)=dq0(ir)+akap*dp0(ir)/r0-akap*p0(ir)/rr2
	    d2q0(ir)=ve*dp0(ir)-akap*dq0(ir)/r0+dve*p0(ir)
     &		                               +akap*q0(ir)/rr2
	    d2p1(ir)=dq1(ir)+akap*dp1(ir)/r0-akap*p1(ir)/rr2
	    d2q1(ir)=ve*dp1(ir)-akap*dq1(ir)/r0-dp0(ir)+dve*p1(ir)
     &		                               +akap*q1(ir)/rr2
	  enddo
	  if(kda.eq.1) then  !! second energy derivatives are needed
	    do ir=1,nrad
	      r0=r(ir)
	      rr2=r0*r0
	      ve=v(ir)-z2/r0-e
	      dve=dv(ir)+z2/rr2
	      dp2(ir)=q2(ir)+akap*p2(ir)/r0
	      dq2(ir)=ve*p2(ir)-akap*q2(ir)/r0-2.d0*p1(ir)
	      d2p2(ir)=dq2(ir)+akap*dp2(ir)/r0-akap*p2(ir)/rr2
	      d2q2(ir)=ve*dp2(ir)-akap*dq2(ir)/r0-2.d0*dp1(ir)+dve*p2(ir)
     &		                                 +akap*q2(ir)/rr2
	    enddo
	  endif
	else if(irel.eq.1) then
	  do ir=1,nrad
	    r0=r(ir)
	    rr2=r0*r0
	    ve=v(ir)-z2/r0-e
          alf=1.d0-ve/c2
	    dve=dv(ir)+z2/rr2
	    dalf=-dve/c2
	    a=akap/alf/rr2
	    dp0(ir)=alf*q0(ir)+p0(ir)/r0
	    dq0(ir)=(ve+a)*p0(ir)-q0(ir)/r0
	    dp1(ir)=alf*q1(ir)+p1(ir)/r0+q0(ir)/c2
	    dq1(ir)=(ve+a)*p1(ir)-q1(ir)/r0-(1.d0+a/alf/c2)*p0(ir)
	    d2p0(ir)=alf*dq0(ir)+dp0(ir)/r0+dalf*q0(ir)-p0(ir)/rr2
	    d2q0(ir)=(ve+a)*dp0(ir)-dq0(ir)/r0
     &		    +(dve-a*(2.d0/r0+dalf/alf))*p0(ir)+q0(ir)/rr2
	    d2p1(ir)=alf*dq1(ir)+dp1(ir)/r0
     &		    +dq0(ir)/c2+dalf*q1(ir)-p1(ir)/rr2
	    d2q1(ir)=(ve+a)*dp1(ir)-dq1(ir)/r0-(1.d0+a/alf/c2)*dp0(ir)
     &		    +(dve-a*(2.d0/r0+dalf/alf))*p1(ir)+q1(ir)/rr2
     &            +2.d0*a/alf/c2*(1.d0/r0+dalf/alf)*p0(ir)
	  enddo
	  if(kda.eq.1) then  !! second energy derivatives are needed
	    do ir=1,nrad
	      r0=r(ir)
	      rr2=r0*r0
	      ve=v(ir)-z2/r0-e
            alf=1.d0-ve/c2
	      dve=dv(ir)+z2/rr2
	      dalf=-dve/c2
	      a=akap/alf/rr2
	      dp2(ir)=alf*q2(ir)+p2(ir)/r0+2.d0*q1(ir)/c2
	      dq2(ir)=(ve+a)*p2(ir)-q2(ir)/r0-2.d0*(1.d0+a/alf/c2)*p1(ir)
     &             +2.d0*a/(alf*c2)**2*p0(ir)
	      d2p2(ir)=alf*dq2(ir)+dp2(ir)/r0+2.d0*dq1(ir)/c2
     &		      +dalf*q2(ir)-p2(ir)/rr2
	      d2q2(ir)=(ve+a)*dp2(ir)-dq2(ir)/r0
     &		      -2.d0*(1.d0+a/alf/c2)*dp1(ir)
     &              +2.d0*a/(alf*c2)**2*dp0(ir)
     &              +(dve-a*(2.d0/r0+dalf/alf))*p2(ir)+q2(ir)/rr2
     &              +4.d0*a/alf/c2*(1.d0/r0+dalf/alf)*p1(ir)
     &              -a/(alf*c2)**2*(4.d0/r0+6.d0*dalf/alf)*p0(ir)
	    enddo
	  endif
	else if(irel.eq.2) then
	  do ir=1,nrad
	    r0=r(ir)
	    rr2=r0*r0
	    ve=v(ir)-z2/r0-e
	    dve=dv(ir)+z2/rr2
          alf=1.d0-ve/c2
	    dalf=-dve/c2
	    dp0(ir)=alf*q0(ir)-akap*p0(ir)/r0
	    dq0(ir)=ve*p0(ir)+akap*q0(ir)/r0
	    dp1(ir)=alf*q1(ir)-akap*p1(ir)/r0+q0(ir)/c2
	    dq1(ir)=ve*p1(ir)+akap*q1(ir)/r0-p0(ir)
	    d2p0(ir)=alf*dq0(ir)-akap*dp0(ir)/r0+akap*p0(ir)/rr2
     &		                                +dalf*q0(ir)
	    d2q0(ir)=ve*dp0(ir)+akap*dq0(ir)/r0+dve*p0(ir)
     &		                               -akap*q0(ir)/rr2
	    d2p1(ir)=alf*dq1(ir)-akap*dp1(ir)/r0+dq0(ir)/c2
     &		    +akap*p1(ir)/rr2+dalf*q1(ir)
	    d2q1(ir)=ve*dp1(ir)+akap*dq1(ir)/r0-dp0(ir)+dve*p1(ir)
     &		                               -akap*q1(ir)/rr2
	  enddo
	  if(kda.eq.1) then  !! second energy derivatives are needed
	    do ir=1,nrad
	      r0=r(ir)
	      rr2=r0*r0
	      ve=v(ir)-z2/r0-e
	      dve=dv(ir)+z2/rr2
            alf=1.d0-ve/c2
	      dalf=-dve/c2
	      dp2(ir)=alf*q2(ir)-akap*p2(ir)/r0+2.d0*q1(ir)/c2
	      dq2(ir)=ve*p2(ir)+akap*q2(ir)/r0-2.d0*p1(ir)
	      d2p2(ir)=alf*dq2(ir)-akap*dp2(ir)/r0+2.d0*dq1(ir)/c2
     &		      +akap*p2(ir)/rr2+dalf*q2(ir)
	      d2q2(ir)=ve*dp2(ir)+akap*dq2(ir)/r0-2.d0*dp1(ir)+dve*p2(ir)
     &		                                 -akap*q2(ir)/rr2
	    enddo
	  endif
	endif
      end
