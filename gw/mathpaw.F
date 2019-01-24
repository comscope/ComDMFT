      subroutine mathpaw(nsitmax,radius,rofi,roatom)
*****************************************************************
*  Perform a summation of the pseudo charge tails over lattice  *
*****************************************************************
c  Note:
c
c    Elements rofi(0,:) and roatom(0,:) are never initialized.
c    Hence the indices in this routine were adjusted to ensure that
c    these undefined elements are never used.
c
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: nsitmax
      real*8, intent(in) :: roatom(0:nrmax,nsort),rofi(0:nrmax,nsort),
     &                           radius
      integer :: isort,nvcmax,ierr,mt,irad,i,i0,isite,jatom,jsort,irad1,
     &               irad2,num,ilmax,ii,l,m,lm,isym,lp,isymp
      real*8 :: sqpi4,rovtt,rmd,rpd,hout1,xpnt
      integer, allocatable :: iref(:,:),nsite(:)
      real*8,allocatable :: work(:,:),answ(:),pol(:),rnew(:),delta(:,:),
     &                          rovt1(:),wd(:,:,:)
      allocate(work(0:nrmax,maxnsym),answ(maxnsym),pol(0:maxpb),
     &             rnew(0:nrmax),delta(nsitmax,nsort),rovt1(0:nrmax),
     &         wd(nsitmax,maxnsym,nsort))
      allocate(iref(nsitmax,nsort),nsite(nsort))
      nvcmax=30000
      sqpi4=2.d0*sqrt(pi)
      answ=0.0d0
      call setsites(delta,iref,nsite,nvcmax,wd,nsitmax,ierr,radius)
      if(ierr.gt.0) then
        if(maswrk) write(iun,*)' Exit from Setsites in Mathpaw : Err'
        call ending
      endif
c *** calculate tail contribution to the pseudo density ***
      do isort=1,nsort
        mt=indmt(1,isort,1)
        call dcopy(nrad(isort)+1,roatom(0,isort),1,ro(mt),1)
        do irad=1,nrad(isort)
          i0=1
          do i=nr(isort),0,-1
            if(rofi(i,isort).lt.r(irad,isort)) then
              i0=max(1,i-1)
              goto 50
            endif
          enddo
50        call fint(roatom(i0,isort),rofi(i0,isort),4,r(irad,isort),
     &              ro(mt+irad))
        enddo
        do isite=1,nsite(isort)
          jatom=iref(isite,isort)
          jsort=is(jatom)
          if(delta(isite,isort).eq.0.d0) goto 80
          i0=1
          do i=1,nr(jsort)
            if(rofi(i,jsort).gt.delta(isite,isort)) then
              i0=max(1,i-2)
              goto 5
            endif
          enddo
5         call fint(roatom(i0,jsort),rofi(i0,jsort),4,
     &              delta(isite,isort),rovtt)
          mt=indmt(1,isort,1)
          ro(mt)=ro(mt)+sqpi4*rovtt*wd(isite,1,isort)
          do irad=1,nrad(isort)
c  make a grid for taking the integrals
            rmd=delta(isite,isort)-r(irad,isort)
            rpd=delta(isite,isort)+r(irad,isort)
            if (rmd.le.0.0d0) then
              if(maswrk) then
                write(iun,*)'ERROR: muffin tin radius is chosen too'
     &                    //' large'
                write(iun,*)'       You chose      :',
     &                      r(nrad(isort),isort)
                write(iun,*)'       Closest contact:',delta(isite,isort)
              endif
              call ending
            endif
            irad1=0
            irad2=1
            do i=1,nr(jsort)
              if(rofi(i,jsort).gt.rmd) then
                irad1=i-1
                goto 6
              endif
            enddo
6           do i=1,nr(jsort)
              if(rofi(i,jsort).gt.rpd) then
                irad2=i
                goto 7
              endif
            enddo
7           num=irad2-irad1+60
c           ilmax=lmpb(isort)
            ilmax=0
            hout1=(rpd-rmd)/num
c  make the integrand
            do i=0,num
              rnew(i)=rmd+i*hout1
              i0=1
              do ii=1,nr(jsort)
                if(rofi(ii,jsort).gt.rnew(i)) then
                  i0=max(1,ii-2)
                  goto 8
                endif
              enddo
8             call fint(roatom(i0,jsort),rofi(i0,jsort),4,rnew(i),
     &                  rovt1(i))
              xpnt=0.5d0*(r(irad,isort)**2+delta(isite,isort)**2-
     &             rnew(i)**2)/(r(irad,isort)*delta(isite,isort))
              call legpol(pol,xpnt,ilmax)
              do l=0,ilmax
                do m=-l,l
                  lm=l*(l+1)+m+1
                  if(sym(lm,isort)) then
                    isym=lmsym(lm,isort)
                    work(i,isym)=rnew(i)*pol(l)
                  endif
                enddo
              enddo
            enddo     !!! over i
            lp=-1
            isymp=0
            do l=0,ilmax
              do m=-l,l
                lm=l*(l+1)+m+1
                if(sym(lm,isort)) then
                  isym=lmsym(lm,isort)
                  if(lp.ne.l) then
                    do i=0,num
                      work(i,isym)=work(i,isym)*rovt1(i)
                    enddo
                    answ(isym)=0.d0
                    do ii=0,num-1
                     answ(isym)=answ(isym)+work(ii,isym)+work(ii+1,isym)
                    enddo
                    answ(isym)=answ(isym)*hout1/2.d0
                    lp=l
                    isymp=isym
                  else
                    answ(isym)=answ(isymp)
                  endif
                endif
              enddo  !!! over m
            enddo     !!! over l
            do l=0,ilmax
              do m=-l,l
                lm=l*(l+1)+m+1
                if(sym(lm,isort)) then
                  isym=lmsym(lm,isort)
                  mt=indmt(isym,isort,1)
                  ro(mt+irad)=ro(mt+irad)+2.d0*pi/r(irad,isort)
     &                       /sqpi4/delta(isite,isort)*answ(isym)*
     &                       wd(isite,isym,isort)
                endif
              enddo
            enddo
          enddo   !!! over irad
80        continue
        enddo   !!! over isite
      enddo   !!! over isort
      deallocate(work,answ,pol,rnew,delta,rovt1,wd,iref,nsite)
      end
