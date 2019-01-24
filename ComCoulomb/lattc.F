      subroutine lattc(dlat0,rlat0,nkdmx,nkrmx,work)
c     sets up the real and reciprocal space lattice vectors
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      integer :: lmax,nkdmx,nkrmx
      double precision, intent(out) :: dlat0(3,nkdmx),rlat0(3,nkrmx)
      double precision :: work(*),tol,tpiba,vol0,
     &  rdist0,qdist0,radd,qadd,a0,tol1,r0,q0
c     ------ set up real and recip vectors ----
      lmax=8
      tol=1.d-65
      tpiba=2.d0*pi/par
      vol0=amega/par**3
      rdist0=vol0**(1.d0/3.d0)
      qdist0=1.d0/rdist0
      radd=.7d0*rdist0
      qadd=.7d0*qdist0
      a0=ewald0/rdist0
      ewald=a0/par
      tol1=tol*par**(lmax+1)
      call lctoff(a0,vol0,lmax,tol1,r0,q0)
c     r0=15.d0
c     q0=9.d0
      call lgen(rbas,r0+radd,nkd,nkdmx,dlat0,work)
      if(maswrk) write(iun,342) r0,r0*par,radd,nkd
 342  format('  r0=',f9.4,'   rc=',f9.4,'   radd=',f9.4,'   nkd=', i7)
      call lgen(gbas,q0+qadd,nkr,nkrmx,rlat0,work)
      if(maswrk) write(iun,341) q0,q0*tpiba,qadd,nkr
 341  format('  q0=',f9.4,'   qc=',f9.4,'   qadd=',f9.4,'   nkr=', i7)
      end
c     -------------------------------------------------------
      subroutine lctoff(a0,v0,lmax,tol,r0,q0)
c     makes limits r0,q0 for sums in real and recip space for a lattice
c     with lattice constant 1.
      use manager_mod
      use parallel_mod
      use units_mod
      implicit none
      integer :: lmax,i
      double precision ::
     $  f(0:10),g(0:10),q1,gq1,a0,v0,tol,r0,q0,gq0,r1,r2,try,q2
      q1=0.001d0
      if(lmax.gt.2) q1=dsqrt(.5d0*(lmax-2))*a0/pi
      gq1=(2d0*pi*q1)**(lmax-2)*dexp(-(pi*q1/a0)**2)*4d0*pi/v0
      if(tol.gt.gq1.and.maswrk) write(iun,*) '**** lctoff: tol gt gq1'
      q2=50.d0
      q0=5.d0
      do i=1,25
        gq0=(2d0*pi*q0)**(lmax-2)*dexp(-(pi*q0/a0)**2)*4d0*pi/v0
        if(gq0.gt.tol) q1=q0
        if(gq0.lt.tol) q2=q0
        q0=.5d0*(q1+q2)
      enddo
c     ---------------------------------------
      r1=0.1d0
      r2=50.d0
      r0=5.d0
      do i=1,25
        call dlmtor(r0,a0,lmax,f,g)
        if(f(lmax).gt.tol) r1=r0
        if(f(lmax).le.tol) r2=r0
        r0=.5d0*(r1+r2)
      enddo
      try=(2d0*pi*q0)**(lmax-2)*dexp(-(pi*q0/a0)**2)*4d0*pi/v0
      end
