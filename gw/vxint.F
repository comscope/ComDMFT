      subroutine vxint(key,vintr0,bvintr0)
************************************************************
*  The calculation of the exchange-correlation potential   *
*  in Interstitial Region.                                 *
************************************************************
      use etot_mod
	use manager_mod
	use solid_mod
	use units_mod
      implicit none
	integer, intent(in) :: key
	complex*16, intent(out) :: bvintr0(3,nplwro)
	complex*16, intent(out) :: vintr0(nplwro,nspin)
      logical :: jimag,jreal,jreal1,jimag1
	integer :: ispin,im,ia,ib,i,j
      real*8 :: gro(3,3),rolap(2),amod(3),gmod(3,3),pi2a,pi4a,const,
     &	        roup,rodn,excup,excdn,vxcup,vxcdn,ddot,rh,absm
	real*8, allocatable :: rhotm(:,:),rhoim(:),rhors(:,:),rhois(:,:),
     &	                     gradro(:,:,:),grad2ro(:,:,:,:),amag(:,:)
	complex*16, allocatable :: absmag(:)
	allocate(rhotm(mmm,nspmag),rhoim(mmm),rhors(mmm,nspmag),
     &	       rhois(mmm,nspmag),gradro(mmm,3,nspmag),
     &         grad2ro(mmm,3,3,nspmag))
	if(magn.eq.2) then
        allocate(amag(mmm,3))
	  allocate(absmag(mmm))
	  absmag=(0.d0,0.d0)
      endif
      pi2a=2.d0*pi/par
      pi4a=pi2a**2
      const=dble(nspmag)/2.d0
      if(magn.eq.1) then
        jimag=.false.
        if(complex_ro) jimag=.true.
        do ispin=1,nspin
          call fourier(rhotm(1,ispin),rhoim,indstar,rointr(1,ispin),
     &                 1,igbsf,gbs,0,0,.true.,jimag)
          if(nloc.gt.0) then
            do ia=1,3
              call fourier(gradro(1,ia,ispin),rhoim,indstar,
     &                     rointr(1,ispin),1,igbsf,gbs,ia,0,.true.,
     &                     jimag)
              do ib=1,3
                call fourier(grad2ro(1,ib,ia,ispin),rhoim,indstar,
     &                       rointr(1,ispin),1,igbsf,gbs,ib,ia,.true.,
     &                       jimag)
              enddo
            enddo
            call dscal(3*mmm,pi2a,gradro(1,1,ispin),1)
            call dscal(9*mmm,pi4a,grad2ro(1,1,1,ispin),1)
          endif
        enddo   !!! over ispin
      else if(magn.eq.2) then
        jreal=.true.
        jimag=.false.
        if(inv.eq.0) jimag=.true.
        call fourier(rhotm,rhoim,indstar,rointr,1,igbsf,gbs,0,0,.true.,
     &               jimag)
        jreal1=.false.
        jimag1=.false.
        if(inv.lt.2) jreal1=.true.
        if(inv.ne.1) jimag1=.true.
        do ia=1,3
          call fourier(amag(1,ia),rhoim,indstar,spintr(ia,:),1,igbsf,
     &                 gbs,0,0,jreal1,jimag1)
        enddo
        do im=1,mmm
          rhors(im,2)=sqrt(ddot(3,amag(im,1),mmm,amag(im,1),mmm))
          rhois(im,2)=0.d0
          rhotm(im,2)=0.5d0*(rhotm(im,1)+rhors(im,2))
          rhotm(im,1)=0.5d0*(rhotm(im,1)-rhors(im,2))
        enddo
        if(nloc.gt.0) then
          call fourier(rhors(1,2),rhois(1,2),indstar,absmag,0,igbsf,gbs,
     &                 0,0,jreal1,jimag1)
          do ia=1,3
            call fourier(gradro(1,ia,1),rhoim,indstar,rointr,1,igbsf,
     &                   gbs,ia,0,.true.,jimag)
            call fourier(gradro(1,ia,2),rhoim,indstar,absmag,1,igbsf,
     &                   gbs,ia,0,jreal1,jimag1)
            do ib=1,3
              call fourier(grad2ro(1,ib,ia,1),rhoim,indstar,rointr,1,
     &                     igbsf,gbs,ib,ia,.true.,jimag)
              call fourier(grad2ro(1,ib,ia,2),rhoim,indstar,absmag,1,
     &                     igbsf,gbs,ib,ia,jreal1,jimag1)
            enddo
          enddo
          call dscal(6*mmm,pi2a,gradro,1)
          call dscal(18*mmm,pi4a,grad2ro,1)
          do ia=1,3
            do im=1,mmm
              rh=0.5d0*(gradro(im,ia,1)+gradro(im,ia,2))
              gradro(im,ia,2)=0.5d0*(gradro(im,ia,1)-gradro(im,ia,2))
              gradro(im,ia,1)=rh
            enddo
            do ib=1,3
              do im=1,mmm
                rh=0.5d0*(grad2ro(im,ib,ia,1)-grad2ro(im,ib,ia,2))
                grad2ro(im,ib,ia,2)=0.5d0*(grad2ro(im,ib,ia,1)+
     &                                     grad2ro(im,ib,ia,2))
                grad2ro(im,ib,ia,1)=rh
              enddo
            enddo
          enddo
        endif
      endif  !!! over magn
      do im=1,mmm
        rodn=rhotm(im,1)*const
        roup=rhotm(im,nspmag)*const
        if(nloc.gt.0) then
          rolap(1)=0.d0
          rolap(2)=0.d0
          do i=1,3
            gro(i,1)=gradro(im,i,1)
            gro(i,2)=gradro(im,i,nspmag)
            gro(i,3)=(gro(i,1)+gro(i,2))*const
            rolap(1)=rolap(1)+grad2ro(im,i,i,1)
            rolap(2)=rolap(2)+grad2ro(im,i,i,nspmag)
          enddo
          amod(1)=sqrt(ddot(3,gro(1,1),1,gro(1,1),1))
          amod(2)=sqrt(ddot(3,gro(1,2),1,gro(1,2),1))
          amod(3)=sqrt(ddot(3,gro(1,3),1,gro(1,3),1))
          do i=1,3
            gmod(i,1)=0.d0
            gmod(i,2)=0.d0
            gmod(i,3)=0.d0
            do j=1,3
              gmod(i,1)=gmod(i,1)+gro(j,1)*grad2ro(im,i,j,1)
              gmod(i,2)=gmod(i,2)+gro(j,2)*grad2ro(im,i,j,nspmag)
              gmod(i,3)=gmod(i,3)+gro(j,3)*(grad2ro(im,i,j,1)+
     &                                      grad2ro(im,i,j,nspmag))
            enddo
            if(amod(1).gt.1.d-10)gmod(i,1)=gmod(i,1)/amod(1)
            if(amod(2).gt.1.d-10)gmod(i,2)=gmod(i,2)/amod(2)
            if(amod(3).gt.1.d-10)gmod(i,3)=gmod(i,3)*const/amod(3)
          enddo
        endif
        call EXCVXC(ROUP,RODN,ilda,0,gro,gmod,amod,rolap,
     &              EXCUP,EXCDN,VXCUP,VXCDN)
        if(key.eq.1) then
          rhors(im,1)=vxcdn
          rhors(im,nspmag)=vxcup
        else
          rhors(im,1)=excdn
          rhors(im,nspmag)=excup
        endif
        rhois(im,1)=0.d0
        rhois(im,nspmag)=0.d0
        if(magn.eq.2) then
          rh=0.5d0*(rhors(im,1)+rhors(im,2))
          rhors(im,2)=0.5d0*(rhors(im,1)-rhors(im,2))
          rhors(im,1)=rh
        endif
      enddo   !!! over im
      if(magn.eq.1) then
        do ispin=1,nspin
          jimag=.false.
          if(complex_ro) jimag=.true.
          call fourier(rhors(1,ispin),rhois(1,ispin),indstar,
     &                 vintr0(1,ispin),0,igbsf,gbs,0,0,.true.,jimag)
        enddo
      else if(magn.eq.2) then
        call fourier(rhors(1,1),rhois(1,1),indstar,vintr0(1,1),0,igbsf,
     &               gbs,0,0,.true.,jimag)
        do im=1,mmm
          absm=sqrt(ddot(3,amag(im,1),mmm,amag(im,1),mmm))
          if(absm.lt.1.d-10) then
            do ia=1,3
              amag(im,ia)=0.d0
            enddo
          else
            do ia=1,3
              amag(im,ia)=amag(im,ia)*rhors(im,2)/absm
            enddo
          endif
        enddo
        jreal=.false.
        jimag=.false.
        if(inv.lt.2) jreal=.true.
        if(inv.ne.1) jimag=.true.
        do ia=1,3
          call dcopy(mmm,0.d0,0,rhois,1)
          call fourier(amag(1,ia),rhois,indstar,bvintr0(ia,:),0,igbsf,
     &                 gbs,0,0,jreal,jimag)
        enddo
      endif   !!! over magn
	do ispin=1,nspin
        call symscal(vintr0(1,ispin))
	enddo
      if(magn.eq.2) then
        jreal=.false.
        jimag=.false.
        if(inv.lt.2) jreal=.true.
        if(inv.ne.1) jimag=.true.
        call symvec(bvintr0,jreal,jimag)
      endif  !!! for magn = 2 only
	deallocate(rhotm,rhoim,rhors,rhois,gradro,grad2ro)	
	if(magn.eq.2) then
        deallocate(amag,absmag)
      endif
      end