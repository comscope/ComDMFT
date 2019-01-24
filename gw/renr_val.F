      subroutine renr_val
****************************************************************
*  Renormalize valence density according to the input charge   *
*  Written by S.Yu.Savrasov (P.N.Lebedev Physical Institute)   *
****************************************************************
	use atom_mod
	use etot_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
      integer :: ispin,isort,mt,irad
      real*8 :: zvtot,anorm,dqdall
	real*8, allocatable :: work1(:)
	allocate(work1(0:maxnrad))
c     &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      zvtot=0.d0
      do ispin=1,nspin
        do isort=1,nsort
          mt=indmt(1,isort,ispin)
          do irad=0,nrad(isort)
            work1(irad)=ro(mt+irad)*weight(irad,isort)/sqrt(4.d0*pi)
          enddo
          zvtot=zvtot
     &         +dqdall(h(isort),work1,nrad(isort))*nhsort(isort)
        enddo   !!! over isort
        zvtot=zvtot+real(rointr(1,ispin))*amegaint
        zvtot=zvtot+amega*dot_product(real(rointr(2:nplwro,ispin)),
     &                                sovr(2:nplwro))
        if(complex_ro) zvtot=zvtot
     &                      -amega
     &                      *dot_product(imag(rointr(2:nplwro,ispin)),
     &	                                 sovi(2:nplwro))
      enddo   !!! over ispin
*  renormalize valence density
      anorm=nelec/zvtot
      if(abs(anorm-1.d0).gt.0.5d0) then
        if(maswrk) write(iun,*)' RENROV : Anorm = ',anorm
c        call ending
      endif
      do ispin=1,nspin
          do isort=1,nsort
            mt=indmt(1,isort,ispin)
	      call dscal(nrad(isort)+1,anorm,ro(mt),1)
          enddo    !!! over isort
      enddo    !!! over ispin
      rointr=anorm*rointr
      if(maswrk) then
        write(iun,1010)nelec
        write(iun,1020)zvtot
        write(iun,1030)anorm
c        print *,' Anorm = ',anorm
      endif
1010  format(1x,'valence charge in whole  elementary cell  must be ',
     &  g14.7)
1020  format(1x,'valence charge  found via one-center expansion is ',
     &  g14.7)
1030  format(1x,'renormalization coefficient of the val.density is ',
     &  g14.7)
	deallocate(work1)
      end