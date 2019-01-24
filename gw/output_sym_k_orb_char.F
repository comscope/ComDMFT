      subroutine output_sym_k_orb_char(ebn,zbn,n0)
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      IMPLICIT none
      integer, intent(in) :: n0
	real*8, intent(in) :: ebn(n0,npnt,nspin)
      complex*16, intent(in) :: zbn(nfun,n0,maxval(ndim3_k),nspin)
      character*3 :: spd(0:4),spdr(9)
      integer :: ispin,k,ibnd,iatom,isort,ndimr,l2m2,li2,jj2,l2,
     &	         mj2,ie2,km2,ie1,km1,ind0,in,jn,in1,jn1,i_len,nn,
     &           i,n,ib,ibn,isp,ind_k
      integer, allocatable :: iref(:),imi(:,:),ima(:,:)
	real*8, allocatable :: pdos0(:,:,:,:,:),d(:),d0(:,:,:)
	data spd/' S ',' P ',' D ',' F ',' G '/
	data spdr/'1/2','1/2','3/2','3/2','5/2','5/2','7/2','7/2','9/2'/
	if(maswrk) then
	  i_len=len_trim(allfile)
	  if(ubi==' gw') then
	    open(3,file=allfile(1:i_len)//'_Orb_Char_Sym_K.gw')
	  else if(ubi==' qp') then
	    open(3,file=allfile(1:i_len)//'_Orb_Char_Sym_K.qp')
	  else if(ubi=='gwg') then
	    open(3,file=allfile(1:i_len)//'_Orb_Char_Sym_K.gwg')
	  else if(ubi=='dft') then
	    open(3,file=allfile(1:i_len)//'_Orb_Char_Sym_K.dft')
	  else if(ubi==' hf') then
	    open(3,file=allfile(1:i_len)//'_Orb_Char_Sym_K.hf')
	  endif
	endif
	nn=min(9,n0)
	allocate(imi(npnt,nspin))
	allocate(ima(npnt,nspin))
	allocate(pdos0(nrel*maxb+1,natom,nn,npnt,nspin))
	allocate(d(n0))
	allocate(d0(n0,npnt,nspin))
	allocate(iref(n0))
      pdos0=0.d0
      imi=0
      ima=0
      d0=0.d0
	do ind_k=1,ndim3_k(me3_k+1)
	  k=n3_mpi_k(me3_k+1)+ind_k
	  do ispin=1,nspin
          isp=min(ispin,nspin_1)
c --- Search for the right bands ------------------------------------
	    n=n_bnd(k,ispin)
	    do ibnd=1,n
	      d0(ibnd,k,ispin)=ebn(ibnd,k,ispin)-chem_pot
	      d(ibnd)=abs(d0(ibnd,k,ispin))
	    enddo
		  call ordering_all(n,d,iref)
c ---------- Check degeneracy --------------------------------
          imi(k,ispin)=10000
          ima(k,ispin)=-10000
          do ib=1,min(n,nn)
            imi(k,ispin)=min(imi(k,ispin),iref(ib))
            ima(k,ispin)=max(ima(k,ispin),iref(ib))
          enddo
c --------- Get PDOS for the given (K;ispin;BND) -----------------
          ibn=0
          do ib=imi(k,ispin),ima(k,ispin)
            ibn=ibn+1
            if(ibn>nn) cycle
            do iatom=1,natom
	        ind0=io_lem(iatom)-1
              isort=is(iatom)
              ndimr=nrel*lmb(isort)+1
              l2m2=0
              do li2=1,ndimr
                if(irel.ne.2) then
		          jj2=2*(li2-1)
	            l2=li2-1
                else if(irel.eq.2) then
		          jj2=li2+mod(li2,2)-1
	            l2=li2/2
	          endif
                do mj2=-jj2,jj2,2
                  l2m2=l2m2+1
                  do ie2=1,ntle(l2,isort)
			        in1=1
			        if(augm(ie2,l2,isort)/='LOC') in1=2
			        do jn1=1,in1
                      km2=ind0+indbasa(jn1,ie2,l2m2,isort)
                      do ie1=1,ntle(l2,isort)
			            in=1
			            if(augm(ie1,l2,isort)/='LOC') in=2
			            do jn=1,in
                          km1=ind0+indbasa(jn,ie1,l2m2,isort)
                          pdos0(li2,iatom,ibn,k,ispin)=
     &                      pdos0(li2,iatom,ibn,k,ispin)
     &					              +conjg(zbn(km1,ib,ind_k,ispin))
     &                                *zbn(km2,ib,ind_k,ispin)
     &                                *ffsmt(jn,jn1,ie1,ie2,li2,
     &                                       isort,isp)
	                  enddo
                      enddo   !!! over ie1
                    enddo
                  enddo   !!! over ie2
                enddo   !!! over mj2
              enddo   !!! over li2
            enddo        !!! over iatom
          enddo               !!! over ib
        enddo
      enddo
	if(nproc_k/=1) then
	  call dgop(pdos0,(nrel*maxb+1)*natom*nn*npnt*nspin,'  +',
     &	          comm_pnt)
	  call igop(imi,npnt*nspin,'  +',comm_pnt)
	  call igop(ima,npnt*nspin,'  +',comm_pnt)
	  call dgop(d0,n0*npnt*nspin,'  +',comm_pnt)
      endif
	do k=1,npnt
	  do ispin=1,nspin
	    if(maswrk) then
	      write(3,*)' '
	      write(3,'(a9,i1,a5,3f8.5,4x,a1)')' Ispin = ',ispin,
     &	                                     ' Q = ',(pnt(i,k),i=1,3),
     & 	                                     sym_label(k)
	      write(3,'(a10,9f7.3)')'ATOM LI  L',
     &	      (d0(ib,k,ispin)*evolt/2.d0,ib=imi(k,ispin),ima(k,ispin))
	      write(3,*)' '
	    endif
	    do iatom=1,natom
            isort=is(iatom)
	      do li2=1,nrel*lmb(isort)+1
              if(irel.ne.2) then
	          l2=li2-1
              else if(irel.eq.2) then
	          l2=li2/2
	        endif
	        if(l2>4) cycle
	        if(maswrk) then
	          if(irel<=1) write(3,'(a4,a6,9f7.4)')txtel(isort),
     &	          spd(l2)//'   ',(pdos0(li2,iatom,ib,k,ispin),ib=1,nn)
	          if(irel==2) write(3,'(a4,a6,9f7.4)')txtel(isort),
     &	          spd(l2)//spdr(li2),(pdos0(li2,iatom,ib,k,ispin),
     &                                ib=1,nn)
              endif
            enddo
          enddo
        enddo                  !!! over ispin
      enddo  !! over k
	deallocate(d,pdos0,iref,imi,ima,d0)
	if(maswrk) close(3)
      END
