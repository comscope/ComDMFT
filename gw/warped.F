      subroutine warped(key,ispin,vint)
c	key=0 - only Hartree part
c	key=1 - Only XC part
c	key=2 - Hartree + XC
	use manager_mod
	use parallel_mod
	use solid_mod
      implicit none
#ifdef MPI
      include 'mpif.h'
#endif
	integer, intent(in) :: key
	complex*16, intent(in) :: vint(nplwro)
	integer :: ispin,iplw,iplw2,ia,ib,ic,ind,ind_plw,i
      integer, allocatable :: n_dim(:),n_mpi(:)
	allocate(n_dim(nproc))
	allocate(n_mpi(nproc))
	call size_shift_par(nplwdif,nproc,n_dim,n_mpi)
	war(:,ispin)=0.d0
      if(complex_ro) wari(:,ispin)=0.d0
      if(key/=0.and.magn.eq.2) then
        if(inv.lt.2) warb=0.d0
        if(inv.ne.1) warbi=0.d0
      endif
      do ind_plw=1,n_dim(me+1)
	  iplw=n_mpi(me+1)+ind_plw
        do iplw2=1,nplwro
          ia=igbs(1,iplw)+igbs(1,iplw2)
          ib=igbs(2,iplw)+igbs(2,iplw2)
          ic=igbs(3,iplw)+igbs(3,iplw2)
          ind=indplw(ia,ib,ic)
          war(iplw,ispin)=war(iplw,ispin)+real(vint(iplw2))*sovr(ind)
          if(complex_ro) then
            war(iplw,ispin)=war(iplw,ispin)-imag(vint(iplw2))*sovi(ind)
            wari(iplw,ispin)=wari(iplw,ispin)
     &                      +imag(vint(iplw2))*sovr(ind)
     &                      +real(vint(iplw2))*sovi(ind)
	    endif
        enddo
        if(key/=0.and.magn.eq.2) then
          do i=1,3
            do iplw2=1,nplwro
              ia=igbs(1,iplw)+igbs(1,iplw2)
              ib=igbs(2,iplw)+igbs(2,iplw2)
              ic=igbs(3,iplw)+igbs(3,iplw2)
              ind=indplw(ia,ib,ic)
              if(inv.lt.2)
     &          warb(i,iplw)=warb(i,iplw)-real(b_intr(i,iplw2))
     &                                   *sovr(ind)
              if(inv.ne.1)
     &          warbi(i,iplw)=warbi(i,iplw)-imag(b_intr(i,iplw2))
     &                                     *sovr(ind)
              if(inv.eq.0) then
                warb(i,iplw)=warb(i,iplw)-imag(b_intr(i,iplw2))
     &                                   *sovi(ind)
                warbi(i,iplw)=warbi(i,iplw)-real(b_intr(i,iplw2))
     &                                     *sovi(ind)
              endif
            enddo
          enddo
        endif
      enddo
      if(nproc/=1) then
#ifdef MPI
        call DGOP(war(1,ispin),nplwdif,'  +',MPI_COMM_WORLD)
        if(complex_ro) call DGOP(wari(1,ispin),nplwdif,'  +',
     &                           MPI_COMM_WORLD)
        if(key/=0.and.magn.eq.2) then
          if(inv.lt.2) call DGOP(warb,3*nplwdif,'  +',MPI_COMM_WORLD)
          if(inv/=1) call DGOP(warbi,3*nplwdif,'  +',MPI_COMM_WORLD)
        endif
#endif        
      endif
	deallocate(n_dim,n_mpi)
      end