      subroutine rcnke0(q,dta,lmax,dl)

      use solid_mod

c     reduced structure constants by ewald method.
c     functions here are real hankel functions without extra factors
c     E=0  q.ne.0

      implicit none
      integer :: lmax
      double precision :: q(3),dta(3)
      complex*16 :: dl((lmax+1)**2)
      call rcnke01(q,dta,lmax,dl)
      call rcnke02(q,dta,lmax,dl)
      end

c     --------------------------------------------
      subroutine rcnke01(q,dta,lmax,dl)
c     k-space part of reduced structure constants (Ewald).
c     tau in alat-units. q and rlat in tpiba-units.

c     ------- We calculate only regular contribution here --------------

      use solid_mod
      use units_mod
      implicit none
      integer :: lmax,nlm,ilm,ir,l,m
      double precision :: q(3),dta(3),r(3),tpi,gamma,fpibv,tpiba,qq,
     &  scalp,r2,r1,den0,sqpi4
      complex*16 :: dl((lmax+1)**2),eiphi,cs
      double precision, allocatable :: yl(:)

      tpi=2.d0*pi
      gamma=0.25d0/(ewald*ewald) !!  = 1/(4*eta^2)
      fpibv=2.d0*tpi/amega      !!  = 4pi/omega
      tpiba=tpi/par             !!  = 2pi/a
      sqpi4=sqrt(4.d0*pi)
      nlm=(lmax+1)**2
      dl=(0.d0,0.d0)
      qq=dot_product(q,q)
      allocate(yl(nlm))
      do ir=1,nkr
        r=tpiba*(q+rlat(:,ir))
        scalp=par*dot_product(r,dta)
        eiphi=dcmplx(dcos(scalp),dsin(scalp))
        r2=dot_product(r,r)
        if(sqrt(r2)>1.d-8) then
          den0=dexp(-gamma*r2)/r2
          r1=dsqrt(r2)
          call sphharm(r(1)/r1,r(2)/r1,r(3)/r1,lmax,yl)
          ilm=0
          do l=0,lmax
            cs=eiphi*den0*r1**l 
c     The above line is = !k+G!^l*e^(i(k+G)delt)*e^(-(k+G)/4eta**2)/!k+G!^2
            do m=-l,l
              ilm=ilm+1
              dl(ilm)=dl(ilm)+yl(ilm)*cs
            enddo
          enddo
        else
          den0=dexp(-gamma*r2)*v_gg_00
          dl(1)=dl(1)+eiphi*den0/sqpi4
        endif
      enddo
      eiphi=(0.d0,-1.d0)*fpibv
      ilm=0
      do l=0,lmax
        eiphi=eiphi*(0.d0,1.d0)
        do m=-l,l
          ilm=ilm+1
          dl(ilm)=eiphi*dl(ilm)
        enddo
      enddo
c     --------- Addings for k=0 (Effective irregular part) ---------------
c     if(qq<1.d-7) then
c     r=par*dta
c     r0=sqrt(dot_product(r,r))
c     ak=(6.d0*pi*pi/amega/dfloat(nqdiv))**(1.d0/3.d0)
c     ------------ L = 0 ----------------------------------------------
c     dl(1)=dl(1)+6.d0*sqrt(pi)/amega/ak**2
c     if(r0<1.d-6) then
c     dl(1)=dl(1)-sqrt(pi)/amega/2.d0/awald**2
c     else if(r0.gt.1.d-6) then
c     ------------ L = 1 ----------------------------------------------
c     call sphharm(r(1)/r0,r(2)/r0,r(3)/r0,lmax,yl)
c     con=4.d0*pi*r0/3.d0/amega
c     do ilm=2,min(nlm,4)
c     dl(ilm)=dl(ilm)-dcmplx(con*yl(ilm),0.d0)
c     enddo
c     endif
c     endif
      deallocate(yl)
      end

c     -------------------------------------------------
      subroutine rcnke02(q,dta,lmax,dl)
c     adds real space part of reduced structure constants (ewald).
      use solid_mod
      use units_mod
      implicit none
      integer :: lmax,ir1,ir,l,m,ilm
      double precision :: r(3),chi(0:40),q(3),dta(3),tpi,srpi,ta,a2,
     &  ta2,cc,r1,r2,gl,qdotr,g0,ddot,erfc_loc
      complex*16 :: cfac,dl((lmax+1)**2),zz
      double precision, allocatable :: yl(:)
      tpi=2.d0*pi
      srpi=dsqrt(pi)
      ta=2.d0*ewald             !! = 2*eta 
      a2=ewald*ewald            !! = eta^2
      ta2=2.d0*a2               !! = 2*eta^2
      cc=4.d0*a2*ewald/srpi     !! = 4*eta^3/sqrt(pi)
      ir1=2
      if(dta(1)**2+dta(2)**2+dta(3)**2.gt.1.d-7) ir1=1
      allocate(yl((lmax+1)**2))
      do ir=ir1,nkd
        r(1)=par*(dlat(1,ir)-dta(1))
        r(2)=par*(dlat(2,ir)-dta(2))
        r(3)=par*(dlat(3,ir)-dta(3))
        r2=ddot(3,r,1,r,1)
        r1=dsqrt(r2)
        call sphharm(r(1)/r1,r(2)/r1,r(3)/r1,lmax,yl)
c     ------ make the xi's from 0 to lmax ------
        chi(0)=erfc_loc(ewald*r1)/r1
        gl=cc*dexp(-a2*r2)/ta2
        do l=1,lmax
          chi(l)=((2*l-1)*chi(l-1)+gl)/r2
          gl=ta2*gl
        enddo 
        qdotr=tpi*(q(1)*dlat(1,ir)+q(2)*dlat(2,ir)+q(3)*dlat(3,ir))
        cfac=dcmplx(dcos(qdotr),dsin(qdotr))
        ilm=0
        do l=0,lmax
          zz=chi(l)*cfac*r1**l
          do m=-l,l
            ilm=ilm+1
            dl(ilm)=dl(ilm)+yl(ilm)*zz
          enddo
        enddo
      enddo
      deallocate(yl)
      if(ir1.eq.1) return
c     --- add dl3 for diagonal structure constants ------
      g0=-ewald/pi
      dl(1)=dl(1)+g0
      end
