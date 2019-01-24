      subroutine exchener(key,ex_en)
c     key = 1 - for integral of V_xc*Ro
c     key = 2 - for integral of E_xc*Ro       
*********************************************************
*  Calculate exchange-correlation contribution to the   *
*            total energy.                              *
*********************************************************
	use atom_mod
	use etot_mod
	use manager_mod
	use solid_mod
	use parallel_mod
      implicit none
#ifdef MPI
      include 'mpif.h'
#endif
	integer, intent(in) :: key
	real*8, intent(out) :: ex_en
      integer :: ishplw,isort,ispin,l,lm,m,isym,mt,irad,i,
     &	         j,ia,ib,ic,ind
      real*8 :: exchint,exchmt,sum,b,dqdall,bxchmt,bxchint,sr,ddot,si
	real*8, allocatable :: work1(:),v1(:),bv1(:)
	complex*16, allocatable :: vintr1(:,:),bvintr1(:,:)
	allocate(work1(0:maxnrad))
      ishplw=nplwro/nproc
      if(ishplw*nproc.lt.nplwro) ishplw=ishplw+1
c  set up full potential to zero, then calculate exchange-correlation
	allocate(v1(maxmt))
      v1=0.d0
	allocate(bv1(maxmtb))
	bv1=0.d0
	allocate(bvintr1(3,nplwro))
	bvintr1=(0.d0,0.d0)
c  calculate density of the exchange-correlation energy
      call vexch(key,v1,bv1)
	allocate(vintr1(nplwro,nspin))
      vintr1=(0.d0,0.d0)
      call vxint(key,vintr1,bvintr1)
c   calculate mt-contribution to the exchange-correlation energy
      exchmt=0.d0
      do isort=1,nsort
        do ispin=1,nspin
          do l=0,lmpb(isort)
            do m=-l,l
              lm=l*(l+1)+m+1
              if(sym(lm,isort))then
                isym=lmsym(lm,isort)
                mt=indmt(isym,isort,ispin)
                do irad=0,nrad(isort)
                  work1(irad)=ro(mt+irad)*v1(mt+irad)*
     &                              dr(irad,isort)*r(irad,isort)**2
                enddo
                sum=dqdall(h(isort),work1,nrad(isort))
                exchmt=exchmt+sum*nhsort(isort)
              endif
            enddo   !!! over m
          enddo   !!! over l
        enddo   !!! over ispin
      enddo   !!! over isort
c  calculate contribution from the interstitial region
      exchint=0.d0
      do ispin=1,nspin
        b=0.d0
        do i=1+ishplw*me,min(nplwro,ishplw*(me+1))
          do j=1,nplwro
            ia=igbs(1,i)+igbs(1,j)
            ib=igbs(2,i)+igbs(2,j)
            ic=igbs(3,i)+igbs(3,j)
            ind=indplw(ia,ib,ic)
            b=b+real(rointr(i,ispin))*real(vintr1(j,ispin))*sovr(ind)
          enddo
        enddo
        if(complex_ro) then
          do i=1+ishplw*me,min(nplwro,ishplw*(me+1))
            do j=1,nplwro
              ia=igbs(1,i)+igbs(1,j)
              ib=igbs(2,i)+igbs(2,j)
              ic=igbs(3,i)+igbs(3,j)
              ind=indplw(ia,ib,ic)
              b=b-imag(rointr(i,ispin))*imag(vintr1(j,ispin))*sovr(ind)
              b=b-(imag(rointr(i,ispin))*real(vintr1(j,ispin))+
     &             real(rointr(i,ispin))*imag(vintr1(j,ispin)))
     &           *sovi(ind)
            enddo
          enddo
        endif
#ifdef MPI
        call DGOP(b,1,'  +',MPI_COMM_WORLD)
#endif        
        b=b*amega
        exchint=exchint+b
      enddo
      exchmt=exchmt
      exchint=exchint
      bxchmt=0.d0
      bxchint=0.d0
      if(magn.eq.2) then
        do isort=1,nsort
          do i=1,3
            do l=0,lmpb(isort)
              do m=-l,l
                lm=l*(l+1)+m+1
                if(symb(lm,i,isort))then
                  isym=lmsymb(lm,i,isort)
                  mt=indmtb(isym,isort)
                  do irad=0,nrad(isort)
                    work1(irad)=spmt(mt+irad)*bv1(mt+irad)
     &                         *dr(irad,isort)*r(irad,isort)**2
                  enddo
                  bxchmt=bxchmt-dqdall(h(isort),work1,nrad(isort))
     & 				           *nhsort(isort)
                endif
              enddo   !!! over m
            enddo   !!! over l
          enddo   !!! over i
        enddo   !!! over isort
c  calculate contribution from the interstitial region
        b=0.d0
        do i=1+ishplw*me,min(nplwro,ishplw*(me+1))
          do j=1,nplwro
            sr=0.d0
            if(inv.lt.2) sr=sr+ddot(3,real(spintr(:,i)),1,
     &                                real(bvintr1(:,j)),1)
            if(inv.ne.1) sr=sr-ddot(3,imag(spintr(:,i)),1,
     &                                imag(bvintr1(:,j)),1)
            ia=igbs(1,i)+igbs(1,j)
            ib=igbs(2,i)+igbs(2,j)
            ic=igbs(3,i)+igbs(3,j)
            ind=indplw(ia,ib,ic)
            b=b-sr*sovr(ind)
            if(inv.eq.0) then
              si=ddot(3,real(spintr(:,i)),1,imag(bvintr1(:,j)),1)+
     &           ddot(3,imag(spintr(:,i)),1,real(bvintr1(:,j)),1)
              b=b+si*sovi(ind)
            endif
          enddo
        enddo
#ifdef MPI
        call DGOP(b,1,'  +',MPI_COMM_WORLD)
#endif        
        b=b*amega
        bxchint=bxchint+b
      endif
      exchmt=exchmt+bxchmt
      exchint=exchint+bxchint
      ex_en=exchmt+exchint
      if(maswrk.and.key==2) then
        write(iun,1010)exchmt,exchint
        if(magn.eq.2) write(iun,1020)bxchmt,bxchint
      endif
1010  format(1x,'(Full): exch-corr.energy : mt,int-parts are',
     &       2(1x,g17.10))
1020  format(1x,'(B):exch-corr.energy : mt,int-parts are',2(1x,g17.10))
	deallocate(work1,v1,vintr1,bv1,bvintr1)
      end