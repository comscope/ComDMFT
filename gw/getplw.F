      subroutine getplw
c     numplw  - number of plane waves in sphere cutoffro+2*(cutoff+akmax)
c     nplwro  - number of plane waves for density (in sphere cutoffro)
c     nplwbas - number of plane waves for basis (in sphere cutoff+akmax)
c     nplwdif - number of plane waves for G'-G in sphere 2*(cutoff+akmax)
	use manager_mod
	use parallel_mod
	use solid_mod
      implicit none
      logical :: err
      integer :: ico,nout1,ncosph,inf,nco_bas,nco_dif,nco_ro,ivec,ig1,
     &           ig2,ig3,i,ig,j,l,k,igb,istar,ibas
      real*8 :: x(3),rout,det,det3,delt,gb
      akmax=0.d0
	do k=0,1
	  do j=0,1
	    do i=0,1
	      x=i*gbas(:,1)+j*gbas(:,2)+k*gbas(:,3)
	      akmax=max(akmax,sqrt(x(1)**2+x(2)**2+x(3)**2))
	    enddo
	  enddo
	enddo
	akmax=akmax/2+1.d-7
      cutoff_g_all=2.d0*max(cutoffro,2*(cutoff+akmax),
     &	                    2*(cutoffgw+akmax))+0.2d0
	allocate(lastg(0:maxplw/2))
	allocate(rcosph(maxplw/2))
	x=0.d0
      call vgensp(cutoff_g_all,maxplw,gbas,x,rout,nout1,gbs,ncosph,
     &            rcosph,lastg(1),inf,err)
      call timel('**** vgensp in getplw finished *****')
      lastg(0)=0
      if(inf.eq.0.and.maswrk) write(iun,1000)nout1,cutoff_g_all
      if(inf.eq.1) then
        if(maswrk) then
          write(iun,1010)rout,nout1
          flush(iun)
        endif
        call ending
      endif
      numplw=nout1
c     &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      det=det3(gbas(1,1),gbas(1,2),gbas(1,3))
      do ivec=1,numplw
         igbs(1,ivec)=idnint(det3(gbs(1,ivec),gbas(1,2),gbas(1,3))/det)
         igbs(2,ivec)=idnint(det3(gbas(1,1),gbs(1,ivec),gbas(1,3))/det)
         igbs(3,ivec)=idnint(det3(gbas(1,1),gbas(1,2),gbs(1,ivec))/det)
      enddo
c -----------------------------------------------------------
      nplwbas=0
      do ico=1,ncosph
         if(rcosph(ico).gt.cutoff+akmax+1.d-13) exit
         nplwbas=lastg(ico)
         nco_bas=ico
      enddo
      nplwdif=0
      do ico=1,ncosph
         if(rcosph(ico).gt.2.d0*(cutoff+akmax+1.d-13)) exit
         nplwdif=lastg(ico)
         nco_dif=ico
      enddo
      nplwdif_pb=0
      do ico=1,ncosph
         if(rcosph(ico).gt.2.d0*(cutoffgw+akmax+1.d-13)) exit
         nplwdif_pb=lastg(ico)
      enddo
      nplwro=0
      do ico=1,ncosph
         if(lastg(ico).gt.maxplwro) then
            if(maswrk) write(iun,1019)maxplwro
            exit
         endif
         if(rcosph(ico).gt.cutoffro+1.d-13) exit
         nplwro=lastg(ico)
         nco_ro=ico
      enddo
      if(maswrk) then
         write(iun,1020)nplwbas
         write(iun,1025)nplwdif
         write(iun,1026)nplwdif_pb
         write(iun,1030)nco_ro
         write(iun,1040)nplwro
      endif
! ------------------------------------------------------------
      maxia=-1000
      maxib=-1000
      maxic=-1000
      do ibas=1,numplw
         maxia=max(maxia,iabs(igbs(1,ibas)))
         maxib=max(maxib,iabs(igbs(2,ibas)))
         maxic=max(maxic,iabs(igbs(3,ibas)))
      enddo
!     &&&&&& Here for density &&&&&&&&&&&&&&&&&&&&&&&
      do ivec=1,nplwro
        ig1=igbs(1,ivec)
        if(ig1.lt.0) ig1=ig1+mdiv(1)
        ig2=igbs(2,ivec)
        if(ig2.lt.0) ig2=ig2+mdiv(2)
        ig3=igbs(3,ivec)
        if(ig3.lt.0) ig3=ig3+mdiv(3)
        igbsf(ivec)=1+ig1+ig2*mdiv(1)+ig3*mdiv(1)*mdiv(2)
      enddo
c     &&&&&&&&&&   Generation of Stars &&&&&&&&&&&&&&&&&&&&&&&&
c     indstar(istar) - last index star=istar
      nstar=0
      indstar(0)=0
      do ico=1,nco_ro
        i=lastg(ico-1)+1
8       nstar=nstar+1
        if(nstar.gt.maxplwro/2) then
          if(maswrk) write(iun,*)' GETPLW: nstar > maxplwro/2'
          call ending
        endif
        rcostar(nstar)=rcosph(ico)
        indstar(nstar)=i
        do ig=1,ngroup
          call rotate(gbs(1,i),gbs(2,i),gbs(3,i),x(1),x(2),x(3),
     &                u(2,ig),2)
          do j=indstar(nstar-1)+1,indstar(nstar)
            delt=abs(gbs(1,j)-x(1))+abs(gbs(2,j)-x(2))+
     &                              abs(gbs(3,j)-x(3))
            if(delt.lt.1.d-5) goto 7
          enddo
          indstar(nstar)=indstar(nstar)+1
          do j=indstar(nstar),lastg(ico)
            delt=abs(gbs(1,j)-x(1))+abs(gbs(2,j)-x(2))+
     &                              abs(gbs(3,j)-x(3))
            if(delt.lt.1.d-5) then
              l=indstar(nstar)
              if(j.ne.l) then
                do k=1,3
                  gb=gbs(k,l)
                  gbs(k,l)=gbs(k,j)
                  gbs(k,j)=gb
                  igb=igbs(k,l)
                  igbs(k,l)=igbs(k,j)
                  igbs(k,j)=igb
                enddo
                igb=igbsf(l)
                igbsf(l)=igbsf(j)
                igbsf(j)=igb
              endif
              goto 7
            endif
          enddo
7         continue
        enddo   !!! over ig
        if(indstar(nstar).lt.lastg(ico)) then
          i=indstar(nstar)+1
          goto 8
        endif
      enddo  !!! over ico
      maxstar=0
      do istar=1,nstar
        maxstar=max(maxstar,indstar(istar)-indstar(istar-1))
      enddo
c ---- Memory distribution for ISTAR-mesh ---------------------------
      call size_shift_par(nstar,nproc,ndim_istar,n_mpi_istar)
      i=n_mpi_istar(me+1) ! zero index for stars
      j=n_mpi_istar(me+1)+ndim_istar(me+1)   ! last index for stars
      nplwrop=indstar(j)-indstar(i)
      if(maswrk) then
        write(iun,1050)nstar
        write(iun,1060)mdiv(1),mdiv(2),mdiv(3)
      endif
      mmm=mdiv(1)*mdiv(2)*mdiv(3)
c ---- Memory distribution for ALL G_Vectors = numplw    ---------
	allocate(ndim_numplw(nproc))
	allocate(n_mpi_numplw(nproc))
      call size_shift_par(numplw,nproc,ndim_numplw,n_mpi_numplw)
      call timel('****** getplw finished *************')
1000  format(1x,'All ',i6,' vectors in sphere ',f9.4,' generated OK')
1010  format(1x,'ERROR from GETPLW : In sphere ',f9.4,' generated ',
     &       i6,' vectors')
1019  format(1x,'WARNING from GETPLW : NUMPLW > MAXPLWRO = ',i9)
1020  format(1x,'Total # of G-vectors for basis found is',i7)
1025  format(1x,'Total # of G"-G  - vectors found is',i7)
1026  format(1x,'Total # of G"-G  for PB is',i7)
1030  format(1x,'Total # of co-spheres for RO and V found is ',i7)
1040  format(1x,'Total # of G-vectors for RO and V found is ',i7)
1050  format(1x,'Total # stars of G-vectors for RO and V found is ',i7)
1060  format(1x,'The divising of the unit cell for density is ',3i5)
1070  format(1x,'The min number of plane waves among proc-s - ',i6,
     &       'max - ',i6)
      end
