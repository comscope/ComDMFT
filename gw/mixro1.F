      subroutine mixro1
*****************************************
*   admixture in the main scf loop      *
*****************************************
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
      integer :: isort,lm,isym,mt,mt1,irad,mti,mti1,i,l,m,istar,i_len
	real*8 :: dz,dm,const,delrho,roout,roinp,ronew,spout,spinp,spnew,
     &	        con,dqdall,ad
      complex*16 :: rooutc,roinpc,ronewc,spoutc,spinpc,spnewc
      real*8, allocatable :: roold(:),work1(:),spmtold(:)
      complex*16, allocatable :: rointro(:,:),spinto(:,:)
      if(iter.eq.0)return
      ad=admix
      if (ubi.ne.'dft') call set_adm(ad)
      allocate(rointro(nplwro,nspin),spinto(3,nplwro))
	allocate(roold(maxmt),work1(0:maxnrad),spmtold(maxmtb))
c open file with the old charge density
      if(maswrk) then
	  i_len=len_trim(allfile)
        open(3,file=allfile(1:i_len)//'_el_density',status='old',
     &       form='unformatted')
        read(3)roold
	  read(3)rointro
      endif
	if(goparr) then
#ifdef MPI
        call brdcst(msgdbl,roold,8*maxmt,master,MPI_COMM_WORLD)
	  call brdcst(msgdbl,rointro,16*nplwro*nspin,master,
     &	            MPI_COMM_WORLD)
#endif     
	endif
      if(magn.eq.2) then
        if(maswrk) then
          read(3)spmtold
	    read(3)spinto
        endif
	  if(goparr) then
#ifdef MPI
          call brdcst(msgdbl,spmtold,8*maxmtb,master,MPI_COMM_WORLD)
	    call brdcst(msgdbl,spinto,48*nplwro,master,MPI_COMM_WORLD)
#endif	    
	  endif
	endif
      if(maswrk) close(3)
      dz=0.d0
      do isort=1,nsort
        isym=1
c define charge convergency between iterations
        mt=indmt(isym,isort,1)
        mt1=indmt(isym,isort,nspin)
        const=dfloat(nspin)/8.d0/pi
        do irad=0,nrad(isort)
          delrho=ro(mt+irad)+ro(mt1+irad)-roold(mt+irad)-
     &                                    roold(mt1+irad)
          work1(irad)=weight(irad,isort)*delrho*delrho
        enddo
        dz=dz+const*dqdall(h(isort),work1(0),nrad(isort))
      enddo   !!! over isort
      if(maswrk) then
        write(iun,1010)dz
        print1010,dz
      endif
      dm=0.d0
      do isort=1,nsort
        do lm=1,(lmpb(isort)+1)**2
          if(sym(lm,isort))then
            isym=lmsym(lm,isort)
c define charge convergency between iterations
            mt=indmt(isym,isort,1)
            mt1=indmt(isym,isort,nspin)
            const=dfloat(nspin)/8.d0/pi
            const=1.d0/4.d0/pi
            do irad=0,nrad(isort)
              delrho=(ro(mt1+irad)-ro(mt+irad))-(roold(mt1+irad)-
     &                                           roold(mt+irad))
              work1(irad)=weight(irad,isort)*delrho*delrho
            enddo
            dm=dm+const*dqdall(h(isort),work1(0),nrad(isort))
c admix new charge density with the old one
            do irad=0,nrad(isort)
              mti=mt+irad
              mti1=mt1+irad
              roout=(ro(mti)+ro(mti1))*dfloat(nspin)/2.d0
              roinp=(roold(mti)+roold(mti1))*dfloat(nspin)/2.d0
              ronew=ad*roout+(1.d0-ad)*roinp
              if(nspin.eq.1) then
                ro(mti)=ronew
              else
                spout=ro(mti1)-ro(mti)
                spinp=roold(mti1)-roold(mti)
                spnew=adspin*spout+(1.d0-adspin)*spinp
                ro(mti1)=0.5d0*(ronew+spnew)
                ro(mti)=0.5d0*(ronew-spnew)
              endif
            enddo   !!! over irad
          endif
        enddo   !!! over lm
        if(magn.eq.2) then
          do i=1,3
            do l=0,lmpb(isort)
              do m=-l,l
                lm=l*(l+1)+m+1
                if(symb(lm,i,isort))then
                  isym=lmsymb(lm,i,isort)
                  mt=indmtb(isym,isort)
                  if(l.eq.0) then
                    do irad=0,nrad(isort)
                      delrho=spmtold(mt+irad)-spmt(mt+irad)
                      work1(irad)=weight(irad,isort)*delrho*delrho
                    enddo
                    dm=dm+dqdall(h(isort),work1(0),nrad(isort))/
     &                            4.d0/pi
                  endif
                  do irad=0,nrad(isort)
                    spmt(mt+irad)=adspin*spmt(mt+irad)+
     &                  (1.d0-adspin)*spmtold(mt+irad)
                  enddo
                endif
              enddo    !!! over m
            enddo    !!! over l
          enddo    !!! over i
        endif  !!! only for magn = 2
      enddo   !!! over isort
      con=dfloat(nspin)/2.d0
      do istar=1,nplwro
        rooutc=(rointr(istar,1)+rointr(istar,nspin))*con
        roinpc=(rointro(istar,1)+rointro(istar,nspin))*con
        ronewc=ad*rooutc+(1.d0-ad)*roinpc
        if(nspin.eq.1) then
          rointr(istar,1)=ronewc
        else
          spoutc=-rointr(istar,1)+rointr(istar,nspin)
          spinpc=-rointro(istar,1)+rointro(istar,nspin)
          spnewc=adspin*spoutc+(1.d0-adspin)*spinpc
          rointr(istar,1)=0.5d0*(ronewc+spnewc)
          rointr(istar,2)=0.5d0*(ronewc-spnewc)
        endif
      enddo
      if(magn.eq.2) then
        do i=1,3
          do istar=1,nplwro
            spintr(i,istar)=spintr(i,istar)*adspin+(1.d0-adspin)*
     &                                      spinto(i,istar)
          enddo
        enddo
      endif  !!! only for magn = 2
      if(maswrk) then
        if(magn.eq.2.or.nspin.eq.2)write(iun,1020)dm
        if(magn.eq.2.or.nspin.eq.2)print 1020,dm
	endif
      deallocate(roold,work1,spmtold)
1010  format(1x,'charge density self-consistency=',g14.7)
1020  format(1x,'magnetization  self-consistency=',g14.7)
      deallocate(rointro,spinto)
      end
