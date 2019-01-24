      subroutine inp_gga(sintet,costet,irad,yl,d12yl,rotot,romod,gro,
     &	               gmod,rolap,drodr,d2rodr2,dspdr,d2spdr2,
     &                   isort,totm,sm)
	use atom_mod
	use manager_mod
	use solid_mod
      implicit none
c	d12yl(1) = d Y / d teta
c	d12yl(2) = d Y / d fi
c	d12yl(3) = d2 Y / d teta2
c	d12yl(4) = d2 Y / d teta / d fi
c	d12yl(5) = d2 Y / d fi2
      integer :: ndimv,isort,irad,ispin,lm,lget,isym,mt,i,isymb,mtb
      real*8 :: yl(limlpb),rotot(2),d12yl(limlpb,5),gro(3,3),
     &	        gmod(3,3),romod(3),rorr(2),ror(2),rot(2),rof(2),rort(2),
     &          rorf(2),rott(2),rotf(2),roff(2),rolap(2),dmdr(3),
     &          dmdt(3),dmdf(3),d2mdr2(3),d2mdt2(3),d2mdf2(3),
     &          d2mdrdt(3),d2mdrdf(3),d2mdtdf(3),totm(3),
     &          drodr(0:nrad(isort),nsym(isort),nspin),
     &          d2rodr2(0:nrad(isort),nsym(isort),nspin),
     &          dspdr(0:nrad(isort),*),
     &          d2spdr2(0:nrad(isort),*),
     &          r1,r2,rst,rst2,rv,sintet,fl,ro0,ro1,ro2,sm,ddot,dsmdr,
     &          dsmdt,dsmdf,d2smdr2,d2smdrdt,d2smdrdf,d2smdt2,d2smdf2,
     &          d2smdtdf,sm_lap,costet,rrr,rrt,rrf,rtt,rtf,rff,rvv,rt,rf
      rv=r(irad,isort)
      r1=1.d0/rv
      r2=r1**2
      rst=r1/sintet
      rst2=rst**2
	ndimv=(lmpb(isort)+1)**2
      ror=0.0d0
      rot=0.0d0
      rof=0.0d0
      rorr=0.0d0
      rort=0.0d0
      rorf=0.0d0
      rott=0.0d0
      rotf=0.0d0
      roff=0.0d0
      rolap=0.0d0
      rotot=0.0d0
	if(nspin.eq.2) then
        do ispin=1,nspin
          ror(ispin)=0.d0
          rot(ispin)=0.d0
          rof(ispin)=0.d0
          rorr(ispin)=0.d0
          rort(ispin)=0.d0
          rorf(ispin)=0.d0
          rott(ispin)=0.d0
          rotf(ispin)=0.d0
          roff(ispin)=0.d0
          rolap(ispin)=0.d0
          rotot(ispin)=0.d0
          do lm=1,ndimv
            if(sym(lm,isort)) then
	        fl=dfloat(lget(lm))
	        fl=fl*(fl+1.d0)
              isym=lmsym(lm,isort)
	        mt=indmt(isym,isort,ispin)
              ro0=ro(mt+irad)
              ro1=drodr(irad,isym,ispin)
              ro2=d2rodr2(irad,isym,ispin)
              rotot(ispin)=rotot(ispin)+ro0*yl(lm)
              rolap(ispin)=rolap(ispin)
     &			        +(ro2+2.d0*r1*ro1-fl*r2*ro0)*yl(lm)
              ror(ispin)=ror(ispin)+ro1*yl(lm)
              rot(ispin)=rot(ispin)+ro0*d12yl(lm,1)
              rof(ispin)=rof(ispin)+ro0*d12yl(lm,2)
              rorr(ispin)=rorr(ispin)+ro2*yl(lm)
              rort(ispin)=rort(ispin)+ro1*d12yl(lm,1)
              rorf(ispin)=rorf(ispin)+ro1*d12yl(lm,2)
              rott(ispin)=rott(ispin)+ro0*d12yl(lm,3)
              rotf(ispin)=rotf(ispin)+ro0*d12yl(lm,4)
              roff(ispin)=roff(ispin)+ro0*d12yl(lm,5)
            endif
          enddo  !!! over lm
        enddo
	else if(nspin.eq.1) then
        ror(1)=0.d0
        rot(1)=0.d0
        rof(1)=0.d0
        rorr(1)=0.d0
        rort(1)=0.d0
        rorf(1)=0.d0
        rott(1)=0.d0
        rotf(1)=0.d0
        roff(1)=0.d0
        rolap(1)=0.d0
        rotot(1)=0.d0
        do lm=1,ndimv
          if(sym(lm,isort)) then
	      fl=dfloat(lget(lm))
	      fl=fl*(fl+1.d0)
            isym=lmsym(lm,isort)
	      mt=indmt(isym,isort,1)
            ro0=ro(mt+irad)/2.d0
            ro1=drodr(irad,isym,1)/2.d0
            ro2=d2rodr2(irad,isym,1)/2.d0
            rotot(1)=rotot(1)+ro0*yl(lm)
            rolap(1)=rolap(1)+(ro2+2.d0*r1*ro1-fl*r2*ro0)*yl(lm)
            ror(1)=ror(1)+ro1*yl(lm)
            rot(1)=rot(1)+ro0*d12yl(lm,1)
            rof(1)=rof(1)+ro0*d12yl(lm,2)
            rorr(1)=rorr(1)+ro2*yl(lm)
            rort(1)=rort(1)+ro1*d12yl(lm,1)
            rorf(1)=rorf(1)+ro1*d12yl(lm,2)
            rott(1)=rott(1)+ro0*d12yl(lm,3)
            rotf(1)=rotf(1)+ro0*d12yl(lm,4)
            roff(1)=roff(1)+ro0*d12yl(lm,5)
          endif
        enddo  !!! over lm
	  if(magn.eq.1) then
c ------------ Copying to ispin=2 ------------------------------
          rotot(2)=rotot(1)
	    rolap(2)=rolap(1)
          ror(2)=ror(1)
	    rot(2)=rot(1)
          rof(2)=rof(1)
	    rorr(2)=rorr(1)
          rort(2)=rort(1)
	    rorf(2)=rorf(1)
          rott(2)=rott(1)
	    rotf(2)=rotf(1)
          roff(2)=roff(1)
c -----------------------------------------------------------
        else if(magn.eq.2) then
c ----------- Calculate derivatives of Magnetization -------------------
          do i=1,3
            dmdr(i)=0.d0
            dmdt(i)=0.d0
            dmdf(i)=0.d0
            d2mdr2(i)=0.d0
            d2mdrdt(i)=0.d0
            d2mdrdf(i)=0.d0
            d2mdtdf(i)=0.d0
            d2mdt2(i)=0.d0
            d2mdf2(i)=0.d0
            totm(i)=0.d0
            do lm=1,ndimv
              if(symb(lm,i,isort)) then
                isymb=lmsymb(lm,i,isort)
	          mtb=indmtb(isymb,isort)
                ro0=spmt(mtb+irad)
                ro1=dspdr(irad,isymb)
                ro2=d2spdr2(irad,isymb)
                totm(i)=totm(i)+ro0*yl(lm)
                dmdr(i)=dmdr(i)+ro1*yl(lm)
                dmdt(i)=dmdt(i)+ro0*d12yl(lm,1)
                dmdf(i)=dmdf(i)+ro0*d12yl(lm,2)
                d2mdr2(i)=d2mdr2(i)+ro2*yl(lm)
                d2mdt2(i)=d2mdt2(i)+ro0*d12yl(lm,3)
                d2mdf2(i)=d2mdf2(i)+ro0*d12yl(lm,5)
                d2mdrdt(i)=d2mdrdt(i)+ro1*d12yl(lm,1)
                d2mdrdf(i)=d2mdrdf(i)+ro1*d12yl(lm,2)
                d2mdtdf(i)=d2mdtdf(i)+ro0*d12yl(lm,4)
              endif
            enddo  !!! over lm
	    enddo
	    sm=sqrt(ddot(3,totm,1,totm,1))
	    if(sm.lt.1.d-12) then
	      dsmdr=0.d0
	      dsmdt=0.d0
	      dsmdf=0.d0
	      d2smdr2=0.d0
	      d2smdrdt=0.d0
	      d2smdrdf=0.d0
	      d2smdt2=0.d0
	      d2smdf2=0.d0
	      d2smdtdf=0.d0
	    else
	      dsmdr=ddot(3,totm,1,dmdr,1)/sm
	      dsmdt=ddot(3,totm,1,dmdt,1)/sm
	      dsmdf=ddot(3,totm,1,dmdf,1)/sm
	      d2smdr2=(ddot(3,dmdr,1,dmdr,1)+ddot(3,totm,1,d2mdr2,1))/sm
     &           -dsmdr**2/sm
	      d2smdrdt=(ddot(3,dmdr,1,dmdt,1)+ddot(3,totm,1,d2mdrdt,1))/sm
     &           -dsmdr*dsmdt/sm
	      d2smdrdf=(ddot(3,dmdr,1,dmdf,1)+ddot(3,totm,1,d2mdrdf,1))/sm
     &           -dsmdr*dsmdf/sm
	      d2smdt2=(ddot(3,dmdt,1,dmdt,1)+ddot(3,totm,1,d2mdt2,1))/sm
     &           -dsmdt**2/sm
	      d2smdf2=(ddot(3,dmdf,1,dmdf,1)+ddot(3,totm,1,d2mdf2,1))/sm
     &           -dsmdf**2/sm
	      d2smdtdf=(ddot(3,dmdt,1,dmdf,1)+ddot(3,totm,1,d2mdtdf,1))/sm
     &           -dsmdt*dsmdf/sm
	    endif
	    sm_lap=2.d0*r1*dsmdr+d2smdr2+costet*r1*rst*dsmdt+r2*d2smdt2
     &              +rst2*d2smdf2
c ------- Forming RO_up and RO_down values ----------------------------
          rotot(2)=rotot(1)-0.5d0*sm
          rolap(2)=rolap(1)-0.5d0*sm_lap
          ror(2)=ror(1)-0.5d0*dsmdr
          rot(2)=rot(1)-0.5d0*dsmdt
          rof(2)=rof(1)-0.5d0*dsmdf
          rorr(2)=rorr(1)-0.5d0*d2smdr2
          rort(2)=rort(1)-0.5d0*d2smdrdt
          rorf(2)=rorf(1)-0.5d0*d2smdrdf
          rott(2)=rott(1)-0.5d0*d2smdt2
          rotf(2)=rotf(1)-0.5d0*d2smdtdf
          roff(2)=roff(1)-0.5d0*d2smdf2
          rotot(1)=rotot(1)+0.5d0*sm
          rolap(1)=rolap(1)+0.5d0*sm_lap
          ror(1)=ror(1)+0.5d0*dsmdr
          rot(1)=rot(1)+0.5d0*dsmdt
          rof(1)=rof(1)+0.5d0*dsmdf
          rorr(1)=rorr(1)+0.5d0*d2smdr2
          rort(1)=rort(1)+0.5d0*d2smdrdt
          rorf(1)=rorf(1)+0.5d0*d2smdrdf
          rott(1)=rott(1)+0.5d0*d2smdt2
          rotf(1)=rotf(1)+0.5d0*d2smdtdf
          roff(1)=roff(1)+0.5d0*d2smdf2
	  endif
	endif
	do ispin=1,2
        gro(1,ispin)=ror(ispin)
        gro(2,ispin)=r1*rot(ispin)
        gro(3,ispin)=rst*rof(ispin)
        romod(ispin)=dot_product(gro(:,ispin),gro(:,ispin))
        romod(ispin)=dsqrt(romod(ispin))
        gmod(1,ispin)=(ror(ispin)*rorr(ispin)+r2*rot(ispin)*
     &      (rort(ispin)-r1*rot(ispin))+rst2*rof(ispin)*
     &      (rorf(ispin)-r1*rof(ispin)))/romod(ispin)
        gmod(2,ispin)=(ror(ispin)*rort(ispin)+r2*rott(ispin)*
     &      rot(ispin)+rst2*rof(ispin)*(rotf(ispin)-
     &      costet/sintet*rof(ispin)))/romod(ispin)
        gmod(3,ispin)=(ror(ispin)*rorf(ispin)+r2*rotf(ispin)
     &		   *rot(ispin)+rst2*rof(ispin)*roff(ispin))/romod(ispin)
        gmod(2,ispin)=r1*gmod(2,ispin)
        gmod(3,ispin)=rst*gmod(3,ispin)
      enddo
      do i=1,3
        gro(i,3)=gro(i,1)+gro(i,2)
      enddo
      romod(3)=0.d0
      do i=1,3
        romod(3)=romod(3)+gro(i,3)**2
      enddo
      romod(3)=dsqrt(romod(3))
      rrr=rorr(1)+rorr(2)
      rrt=rort(1)+rort(2)
      rrf=rorf(1)+rorf(2)
      rtt=rott(1)+rott(2)
      rtf=rotf(1)+rotf(2)
      rff=roff(1)+roff(2)
      rvv=ror(1)+ror(2)
      rt=rot(1)+rot(2)
      rf=rof(1)+rof(2)
      gmod(1,3)=(rvv*rrr+r2*rt*(rrt-r1*rt)+rst2*rf*(rrf-r1*rf))
     &	                                      /romod(3)
      gmod(2,3)=(rvv*rrt+r2*rtt*rt+rst2*rf*(rtf-costet/sintet*rf))
     &                                                    /romod(3)
      gmod(3,3)=(rvv*rrf+r2*rtf*rt+rst2*rf*rff)/romod(3)
      gmod(2,3)=r1*gmod(2,3)
      gmod(3,3)=rst*gmod(3,3)
      end
