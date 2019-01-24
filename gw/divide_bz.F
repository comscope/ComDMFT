      subroutine divide_bz
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer :: ia,ib,ic,m,ierror
	real*8 :: x1,y1,z1,rad_inscribed_sphere,qb(3,3),ct0
2	nr_full=nrdiv(1)*nrdiv(2)*nrdiv(3)
      if (nr_full<=0) then
        if (maswrk) then
          write(iun,*)'divide_bz: Something bad happened to the grid'
          write(iun,*)'divide_bz: nrdiv = ',nrdiv
        endif
        call ending
        stop
      endif
      r0b0(:,1) = rbas(:,1)/nrdiv(1)
      r0b0(:,2) = rbas(:,2)/nrdiv(2)
      r0b0(:,3) = rbas(:,3)/nrdiv(3)
      q0b0(:,1) = gbas(:,1)*nrdiv(1)
      q0b0(:,2) = gbas(:,2)*nrdiv(2)
      q0b0(:,3) = gbas(:,3)*nrdiv(3)
	if(allocated(r_pnt)) deallocate(r_pnt)
	if(allocated(r_pnt_ind)) deallocate(r_pnt_ind)
      allocate(r_pnt_ind(nr_full))
      allocate(r_pnt(3,nr_full))
c ------------ Here we set up r-mesh for given NRDIV     ----------
	call get_r(ierror)
	if(ierror/=0) goto 2
c --------- Setting up CUTOFF and CUTOFFGW for given NRDIV  ---
      qb(:,1) = gbas(:,1)*mdiv(1)
      qb(:,2) = gbas(:,2)*mdiv(2)
      qb(:,3) = gbas(:,3)*mdiv(3)
	cutoffro=rad_inscribed_sphere(qb)
	ct0=rad_inscribed_sphere(q0b0)
	cutoffgw=cut_pb_ratio*ct0
	cutoff=cut_lapw_ratio*ct0
	if(maswrk) then
	  write(iun,'(a10,f8.4,a12,f8.4,a12,f8.4)')' CUTOFF = ',cutoff,
     &	               ' CUTOFFRO = ',cutoffro,' CUTOFFGW = ',cutoffgw
      endif
c --------------------------------------------------------------
12	nqdiv=ndiv(1)*ndiv(2)*ndiv(3)
      nttrmaxs=nqdiv*6
	if(allocated(pnt)) deallocate(pnt)
	if(allocated(wgt)) deallocate(wgt)
	if(allocated(iwtet)) deallocate(iwtet)
	if(allocated(idtet)) deallocate(idtet)
      allocate(pnt(3,nqdiv))
      allocate(wgt(nqdiv))
      allocate(iwtet(nttrmaxs))
      allocate(idtet(4,nttrmaxs))
	do m = 1, 3
        rb0(m,1) = rbas(m,1)*ndiv(1)
        rb0(m,2) = rbas(m,2)*ndiv(2)
        rb0(m,3) = rbas(m,3)*ndiv(3)
        qb0(m,1) = gbas(m,1)/ndiv(1)
        qb0(m,2) = gbas(m,2)/ndiv(2)
        qb0(m,3) = gbas(m,3)/ndiv(3)
      enddo
      call get_k(ierror)
	if(ierror/=0) then
	  ndiv(ierror)=ndiv(ierror)+1
	  if(maswrk) write(iun,'(a34,3i5)')
     &	            ' Error in GET_K, NDIV is set to be',
     &                            (ndiv(m),m=1,3)
	  goto 12
	endif
c ---------------------------------------------------------------------
	if(maswrk) then
	  write(iun,'(a7,3i4)')' NDIV: ',(ndiv(m),m=1,3)
	  write(iun,'(a7,3i4)')'NRDIV: ',(nrdiv(m),m=1,3)
	  write(iun,'(a9,i4,a9,i4)')'NR_FULL: ',nr_full,' NR_RED: ',nr_red
	endif
	do m=1,3
	  ndiv_full(m)=ndiv(m)*nrdiv(m)
	enddo
	r_per(:,1)=rbas(:,1)*ndiv(1)
	r_per(:,2)=rbas(:,2)*ndiv(2)
	r_per(:,3)=rbas(:,3)*ndiv(3)
	g_per(:,1)=gbas(:,1)*nrdiv(1)
	g_per(:,2)=gbas(:,2)*nrdiv(2)
	g_per(:,3)=gbas(:,3)*nrdiv(3)
	cut_rr_eff=rad_inscribed_sphere(r_per)
c --------- Getting the minimal distances between vectors of lattices --
	rbas_min=1000.d0
	do ic=-ndiv(3),ndiv(3),ndiv(3)
	  do ib=-ndiv(2),ndiv(2),ndiv(2)
	    do ia=-ndiv(1),ndiv(1),ndiv(1)
	      if(ia==0.and.ib==0.and.ic==0) cycle
	      x1=rbas(1,1)*ia+rbas(1,2)*ib+rbas(1,3)*ic
	      y1=rbas(2,1)*ia+rbas(2,2)*ib+rbas(2,3)*ic
	      z1=rbas(3,1)*ia+rbas(3,2)*ib+rbas(3,3)*ic
	      rbas_min=min(x1*x1+y1*y1+z1*z1,rbas_min)
	    enddo
	  enddo
	enddo
	rbas_min=sqrt(rbas_min)
	gbas_min=1000.d0
	do ic=-nrdiv(3),nrdiv(3),nrdiv(3)
	  do ib=-nrdiv(2),nrdiv(2),nrdiv(2)
	    do ia=-nrdiv(1),nrdiv(1),nrdiv(1)
	      if(ia==0.and.ib==0.and.ic==0) cycle
	      x1=gbas(1,1)*ia+gbas(1,2)*ib+gbas(1,3)*ic
	      y1=gbas(2,1)*ia+gbas(2,2)*ib+gbas(2,3)*ic
	      z1=gbas(3,1)*ia+gbas(3,2)*ib+gbas(3,3)*ic
	      gbas_min=min(x1*x1+y1*y1+z1*z1,gbas_min)
	    enddo
	  enddo
	enddo
	gbas_min=sqrt(gbas_min)
      call timel('*** Divide_BZ   finished ***********')
      end
