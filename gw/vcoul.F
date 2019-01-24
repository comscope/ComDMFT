      subroutine vcoul
************************************************************
*  The calculation of the coulomb potential within MTs     *
*  without nuclear part.                                   *
*  Written by S.Yu.Savrasov (P.N.Lebedev Physical Instit.) *
************************************************************
	use atom_mod
	use etot_mod
      use manager_mod
      use parallel_mod
	use solid_mod
	use units_mod
      implicit none
      integer :: isort,l,m,lm,isym,mt,mt1,irad
      real*8 :: sqpi4,sl,pl,ctl
	real*8, allocatable :: vc1(:),vc2(:),tmp(:),rotot(:)
	allocate(vc1(0:maxnrad),vc2(0:maxnrad))
	allocate(tmp(0:maxnrad),rotot(0:maxnrad))
c  **** calculate coulomb contribution to the potential ****
      sqpi4=sqrt(4.d0*pi)
	if(iheg) then
	  v_mt_h=0.d0
	  v_intr_h=(0.d0,0.d0)
	else
        do isort=1,nsort
c  do over spherical harmonics expansions
          do l=0,lmpb(isort)
            sl=smt(isort)**l
            pl=2.d0/(2.d0*l+1.d0)
            do m=-l,l
              lm=l*(l+1)+m+1
              if(sym(lm,isort))then
                isym=lmsym(lm,isort)
                mt=indmt(isym,isort,1)
                mt1=indmt(isym,isort,nspin)
                v_mt_h(mt:mt+nrad(isort))=0.d0
                do irad=0,nrad(isort)
                  rotot(irad)=ro(mt+irad)
                  if(nspin==2) rotot(irad)=rotot(irad)+ro(mt1+irad)
                enddo
c   coulomb contribution  (internal part)
                call poiss1(vc1(0),vc2(0),l,nrad(isort),h(isort),
     &			            r(0,isort),weight(0,isort),rl(0,l,isort),
     &                      rotot)
                do irad=0,nrad(isort)
                  v_mt_h(mt+irad)=v_mt_h(mt+irad)
     &				             +pl*(vc1(irad)+vc2(irad))
                enddo
c   coulomb contribution  (external part)
                ctl=vatmt(isym,isort)-pl*vc1(nrad(isort))
                if(l.eq.0) ctl=ctl+2.d0*sqpi4*z(isort)/smt(isort)
                ctl=ctl/sl
                do irad=0,nrad(isort)
                  v_mt_h(mt+irad)=v_mt_h(mt+irad)+ctl*rl(irad,l,isort)
                enddo
              endif
            enddo    !!! over m
          enddo    !!! over l
        enddo    !!! over isort
	endif
      call zerov
c -------- Calculation the Z_VNUCL ---------------------------------------
      z_vnucl=0.d0
      do isort=1,nsort
        mt=indmt(1,isort,1)
        z_vnucl=z_vnucl-0.5d0*z(isort)*v_mt_h(mt)/sqpi4
     &	                                       *dfloat(nhsort(isort))
	enddo
	deallocate(vc1,vc2,tmp,rotot)
      call timel('****   vcoul finished **************')
      end