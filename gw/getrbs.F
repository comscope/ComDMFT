      subroutine getrbs
	use manager_mod
	use parallel_mod
	use solid_mod
      implicit none
      logical :: err
      integer :: nout1,ncosph,inf,ivec,ic,ib,ia,nmax
      real*8 :: rout,det,det3,x1,y1,z1,v(3)
	integer, allocatable :: lst(:)
	real*8, allocatable :: rcs(:),rbs0(:,:)
c --- The linear size of active polygon ----	
	cut_rr=0.d0
	do ic=-ndiv(3),ndiv(3),ndiv(3)
	  do ib=-ndiv(2),ndiv(2),ndiv(2)
	    do ia=-ndiv(1),ndiv(1),ndiv(1)
	      x1=rbas(1,1)*ia+rbas(1,2)*ib+rbas(1,3)*ic
	      y1=rbas(2,1)*ia+rbas(2,2)*ib+rbas(2,3)*ic
	      z1=rbas(3,1)*ia+rbas(3,2)*ib+rbas(3,3)*ic
	      cut_rr=max(x1*x1+y1*y1+z1*z1,cut_rr)
	    enddo
	  enddo
	enddo
	cut_rr=sqrt(cut_rr)
	nmax=maxplw
1     allocate(rbs0(3,nmax))
	allocate(lst(0:nmax/2))
	allocate(rcs(nmax/2))
	v=0.d0
      call vgensp(cut_rr+1.d-6,nmax,rbas,v,rout,nout1,rbs0,ncosph,rcs,
     &            lst(1),inf,err)
      if(err) then
        deallocate(rbs0,lst,rcs)
        nmax=2*nmax
        goto 1
      endif
	deallocate(lst,rcs)
      call timel('**** vgensp in getrbs finished *****')
      if(inf.eq.0.and.maswrk) write(iun,1000)nout1,cut_rr+1.d-6
      if(inf.eq.1) then
        if(maswrk) write(iun,1010)rout,nout1
        call ending
      endif
      n_rr_max=nout1
c     &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      allocate(rbs(3,n_rr_max))
      rbs=rbs0(:,1:n_rr_max)
      deallocate(rbs0)
      allocate(irbs(3,n_rr_max))
      det=det3(rbas(1,1),rbas(1,2),rbas(1,3))
      do ivec=1,n_rr_max
         irbs(1,ivec)=idnint(det3(rbs(1,ivec),rbas(1,2),rbas(1,3))/det)
         irbs(2,ivec)=idnint(det3(rbas(1,1),rbs(1,ivec),rbas(1,3))/det)
         irbs(3,ivec)=idnint(det3(rbas(1,1),rbas(1,2),rbs(1,ivec))/det)
      enddo
	maxra=-1000
      maxrb=-1000
      maxrc=-1000
      do ivec=1,n_rr_max
         maxra=max(maxra,iabs(irbs(1,ivec)))
         maxrb=max(maxrb,iabs(irbs(2,ivec)))
         maxrc=max(maxrc,iabs(irbs(3,ivec)))
      enddo
c --------------------------------------------------------------
      allocate(ind_rsp(-maxra:maxra,-maxrb:maxrb,-maxrc:maxrc))
	ind_rsp=0
	do ivec=1,n_rr_max
        ind_rsp(irbs(1,ivec),irbs(2,ivec),irbs(3,ivec))=ivec
      enddo
c --------------------------------------------------------------
      call timel('****** getrbs finished *************')
c     &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
1000  format(1x,' GETRBS: All ',i9,' vectors in sphere ',f9.4,
     &       ' generated OK')
1010  format(1x,'WARNING from GETPLW : In sphere ',f9.4,' generated ',
     &       i9,' vectors')
      end