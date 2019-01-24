      subroutine vexch(key,v0,bv0)
************************************************************
*  The calculation of the exchange-correlation potential   *
*  within MTs using direct spherical-harmonic expansion.   *
*  Adapted by S.Savrasov from Methfessel program.          *
*  Adapted by A.Kutepov from S.Savrasov program.           *
*  key=1 - potential, =2 - energy density                  *
*  last revision 14.10.06.                                 *
************************************************************
	use atom_mod
      use etot_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
#ifdef MPI
      include 'mpif.h'
#endif
	integer, intent(in) :: key
	real*8, intent(out) :: v0(maxmt),bv0(maxmtb)
	integer :: isort,ndimv,ntheta,nphi,nangl,ishang,iang,itheta,iphi,
     &	         irad,ispin,l,m,lm,isym,mt,i
	real*8 :: dphi,zcoord,xcoord,ycoord,stheta,wgtl,phi,excup,excdn,
     &	        vxcup,vxcdn,wv,sm,rotot(2),vxcpot,bpot,bx,by,bz,
     &          gro(3,3),gmod(3,3),romod(3),rolap(2),vx(2),
     &	        totm(3)
	real*8, allocatable :: zl(:),wl(:),vt(:),bt(:)
	real*8, allocatable :: ylm(:)
	allocate(zl(4*maxpb+3),wl(4*maxpb+3),vt(maxmt))
	vt=0.d0
	if(magn==2) then
	  allocate(bt(maxmtb))
	  bt=0.d0
	endif
	allocate (ylm(limlpb))
c  **** calculate exchange-correlation by direct expansion ****
      do isort=1,nsort
        ndimv=(lmpb(isort)+1)**2
        ntheta=4*lmpb(isort)+3
        nphi=ntheta
        nangl=nphi*ntheta
        ishang=nangl/nproc
        if(ishang*nproc.lt.nangl) ishang=ishang+1
        call mklegw(ntheta,zl,wl)
        dphi=2.d0*pi/dble(nphi)
        do iang=1+ishang*me,min(nangl,ishang*(me+1))
          itheta=iang/nphi
          if(mod(iang,nphi).ne.0) itheta=itheta+1
          iphi=iang-nphi*(itheta-1)
          zcoord=zl(itheta)
          stheta=sqrt(1.d0-zcoord*zcoord)
          wgtl=dphi*wl(itheta)
          phi=dble(iphi-1)*dphi
          xcoord=cos(phi)*stheta
          ycoord=sin(phi)*stheta
          call sphharm(xcoord,ycoord,zcoord,lmpb(isort),ylm)
          do irad=1,nrad(isort)
c  Set up total density at the point "r,theta,phi"
			call inp_lda(irad,ylm,rotot,isort,totm,sm)
c  find contribution to exchange-correlation potential from "r,theta,phi"
            call EXCVXC(rotot(2),rotot(1),ilda,0,gro,gmod,romod,
     &		            rolap,EXCUP,EXCDN,VXCUP,VXCDN)
            if(key.eq.1) then
              vx(1)=vxcdn
              vx(2)=vxcup
            else
              vx(1)=excdn
              vx(2)=excup
            endif
            if(magn==1) then
              do ispin=1,nspin
                wv=wgtl*vx(ispin)
                do l=0,lmpb(isort)
                  do m=-l,l
                    lm=l*(l+1)+m+1
                    if(sym(lm,isort))then
                      isym=lmsym(lm,isort)
                      mt=indmt(isym,isort,ispin)
                      vt(mt+irad)=vt(mt+irad)+ylm(lm)*wv
                    endif
                  enddo   !!! over m
                enddo   !!! over l
              enddo   !!! over ispin
            else if(magn.eq.2) then
              vxcpot=0.5d0*(vx(1)+vx(2))
              wv=wgtl*vxcpot
              do l=0,lmpb(isort)
                do m=-l,l
                  lm=l*(l+1)+m+1
                  if(sym(lm,isort))then
                    isym=lmsym(lm,isort)
                    mt=indmt(isym,isort,1)
                    vt(mt+irad)=vt(mt+irad)+ylm(lm)*wv
                  endif
                enddo   !!! over m
              enddo   !!! over l
              bpot=0.5d0*(vx(1)-vx(2))
              if(sm.lt.1.d-10) then
                bx=0.d0
                by=0.d0
                bz=0.d0
              else
                bx=bpot*totm(1)/sm
                by=bpot*totm(2)/sm
                bz=bpot*totm(3)/sm
              endif
              do i=1,3
                if(i.eq.1) wv=wgtl*bx
                if(i.eq.2) wv=wgtl*by
                if(i.eq.3) wv=wgtl*bz
                do l=0,lmpb(isort)
                  do m=-l,l
                    lm=l*(l+1)+m+1
                    if(symb(lm,i,isort))then
                      isym=lmsymb(lm,i,isort)
                      mt=indmtb(isym,isort)
                      bt(mt+irad)=bt(mt+irad)+wv*ylm(lm)
                    endif
                  enddo   !!! over m
                enddo   !!! over l
              enddo   !!! over i
            endif  !!! only for mag = 2
          enddo   !!! over irad
        enddo   !!! over iang
      enddo   !!! over isort
#ifdef MPI
      call DGOP(vt,maxmt,'  +',MPI_COMM_WORLD)
      if(magn.eq.2) call DGOP(bt,maxmtb,'  +',MPI_COMM_WORLD)
#endif      
      do i=1,maxmt
        v0(i)=v0(i)+vt(i)
      enddo
      if(magn.eq.2) then
        do i=1,maxmtb
          bv0(i)=bv0(i)+bt(i)
        enddo
      endif
      do isort=1,nsort
        do ispin=1,nspin
          mt=indmt(1,isort,ispin)
          v0(mt)=v0(mt+1)
        enddo
        if(magn.eq.2) then
          do i=1,3
            if(symb(1,i,isort)) then
              isym=lmsymb(1,i,isort)
              mt=indmtb(isym,isort)
              bv0(mt)=bv0(mt+1)
            endif
          enddo
        endif
      enddo   !!! over isort
	deallocate (ylm)
	deallocate(zl,wl,vt)
	if(magn==2) deallocate(bt)
      end