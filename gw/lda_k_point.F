      subroutine lda_k_point(ind_k,ispin,kp,err_val,err_slo,ebn,nbn,zbn,
     &                       evbn,indg,nbas)
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer, intent(in) :: ispin,ind_k,indg(nbndtm),nbas
	real*8, intent(in) :: kp(3)
	real*8, intent(inout) :: err_val,err_slo
	integer, intent(out) :: nbn
	real*8, intent(out) :: ebn(nbndf)
	complex*16, intent(out) :: zbn(nfun,nbndf),evbn(nbasmpw,nbndf)
	integer :: nev,neig,ibnd,ind1,iatom,isort,nlct,ndimb,i0,ie,ndim,j,
     &	         i,iv,k
	real*8 :: eps_degen,vabs,v(3),pi2a
      integer, allocatable :: ind_y(:,:)
	real*8, allocatable :: yr(:,:),yi(:,:),ar0(:,:),diag(:),diagb(:),
     &	                     ar1(:,:),br0(:,:),br1(:,:),ai0(:,:),
     &                       ai1(:,:),bi0(:,:),bi1(:,:),gb(:,:),ac(:,:)
	complex*16, allocatable :: a(:,:),b(:,:),yc(:,:),ztmp(:,:),
     &	                         z_tmp(:),sig_x0(:,:)
      k=n3_mpi_k(me3_k+1)+ind_k
      eps_degen=1.d-7
	allocate(gb(3,nplwbas))
      do i=1,nbas/nrel
	  gb(:,i)=gbs(:,indg(i))
	enddo
	ndim=nbas
	pi2a=(pi+pi)/par
      ndim=nbas+nlcloc
      allocate(diag(ndim),diagb(ndim))
      allocate(a(ndim,ndim),b(ndim,ndim))
      allocate(yc(max_bs,nbas))
	a=(0.d0,0.d0)
	b=(0.d0,0.d0)
      allocate(yr(max_bs,nbas))
      allocate(yi(max_bs,nbas))
	yr=0.d0
	yi=0.d0
      allocate(ind_y(2*limlb,natom))
      call setk2(kp,gb,yr,yi,nbas,ind_y)
	yc=dcmplx(yr,yi)
      allocate(ar1(nbas,nbas),br1(nbas,nbas),ai1(nbas,nbas),
     &         bi1(nbas,nbas))
      call setwmt(gb,kp,nbas,ar1,ai1,br1,bi1,indg,ispin,1)
      allocate(ar0(nbas,nlcloc),br0(nbas,nlcloc),ai0(nbas,nlcloc),
     &         bi0(nbas,nlcloc))
c     &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c      call setdiag(gb,b,a,kp,nbas,ndim,yr,yi,ind_y,ar0,ai0,br0,bi0,
c     &             ar1,ai1,br1,bi1,ispin,indg)
      call set_h0(b,a,nbas,ndim,yr,yi,ind_y,ar0,ai0,br0,bi0,
     &            ar1,ai1,br1,bi1,ispin)
      call eiglin(ispin,nbas,a,b,ndim,yr,yi,ind_y,ar0,ai0,br0,bi0,
     &            ar1,ai1,br1,bi1)
      deallocate(ar0,ar1,ai0,ai1,br0,br1,bi0,bi1)
	do j=1,ndim
	  do i=1,j-1
	    a(i,j)=conjg(a(j,i))
	    b(i,j)=conjg(b(j,i))
	  enddo
	enddo
	if(ubi/='dft') then
	  allocate(sig_x0(ndim,ndim))
	  call sigma_x_loc(nbas,ndim,sig_x0,yc,ind_y)
	  call sigma_x_mt_mt(nbas,ndim,sig_x0,yc,ind_y,pvv(1,1,ind_k))
        call sigma_x_mt_int(nbas,ndim,sig_x0,yc,ind_y,pvg(1,1,ind_k))
	  call sigma_x_int_int(nbas,ndim,sig_x0,pgg(1,1,ind_k))
	  a=a+sig_x0
	endif
	allocate(ztmp(ndim,ndim))
c -------- Spectrum of complex matrices ---------------------------
	if(me_tau_kk_pbr==0) call eig_val_gener_solver(eps_degen,ndim,
     & 	                                             ndim,a,b,nev,diag,
     &                                               diagb)
      if(nproc_tau/=1) then
	  call brdcst(msgdbl,a,16*ndim**2,0,comm_tau_kk_pbr)
	  call brdcst(msgdbl,diag,8*ndim,0,comm_tau_kk_pbr)
	  call brdcst(msgdbl,diagb,8*ndim,0,comm_tau_kk_pbr)
	  call brdcst(msgint,nev,4,0,comm_tau_kk_pbr)
      endif
	ztmp=a
	neig=nev
	if(maswrk.and.iprn>=1) then
	  write(iun,*)'ndim new e_overl(1) ',ndim,nev,diagb(1)
	endif
	deallocate(a,b)
c     &&&&&& Getting New N_BND(k) &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	if(nbndf>=neig) then
	  nbn=neig
	else if(nbndf<neig) then
	  do ibnd=nbndf,1,-1
	    if(diag(ibnd+1)-diag(ibnd).gt.1.d-5)then
	      nbn=ibnd
	      exit
	    endif
	  enddo
      endif
	neig=nbn
	ebn(1:nbn)=diag(1:nbn)
      if(nbndf>nbn) ebn(nbn+1:nbndf)=ebn(nbn)+0.2d0
	deallocate(diag,diagb)
	if(ubi/='dft') then
        call sigma_x_final(ndim,sig_x0,ztmp,nbn,ind_k,
     &		               ispin)
	  deallocate(sig_x0)
	endif
c -- Copying first NLCLOC components of EV_LDA into Z for LO - basis --
      call zlacpy(' ',nlcloc,neig,ztmp,ndim,zbn,nfun)
c  ------- Forming rest components of Z from PW-part of basis ---------
      call zgemm('n','n',max_bs,neig,nbas,(1.d0,0.d0),yc,max_bs,
     &           ztmp(nlcloc+1,1),ndim,(0.d0,0.d0),zbn(nlcloc+1,1),nfun)
	deallocate(yc)
c - Transformation ZC from FLAPW representation to MB representation -
      allocate(z_tmp(nfun))
      do ibnd=1,neig
        z_tmp=zbn(:,ibnd)
c --------- Local orbitals -------------------------------------
        ind1=1
        do iatom=1,natom
          isort=is(iatom)
          nlct=nlctot(isort)
	    call zcopy(nlct,z_tmp(ind1),1,zbn(io_lem(iatom),ibnd),1)
	    ind1=ind1+nlct
	  enddo
c --------- Augmented orbitals -------------------------------------
        ind1=nlcloc+1
        do iatom=1,natom
          isort=is(iatom)
          ndimb=2*nrel*(lmb(isort)+1)**2
	    nlct=nlctot(isort)
	    i0=io_lem(iatom)+nlct
	    call zcopy(ndimb,z_tmp(ind1),1,zbn(i0,ibnd),1)
	    ind1=ind1+ndimb
	  enddo
	enddo
	deallocate(z_tmp)
c     &&&&& Forming EV_LDA &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      call zlacpy(' ',nbas,neig,ztmp(nlcloc+1,1),ndim,evbn,nbasmpw) 
	deallocate(ztmp)
c ------------------------------------------------------------------
	allocate(ac(3,nbn))
	ac=0.d0
	if(irel/=2.or..not.rel_interst) 
     &	call check_band_states(zbn,evbn,kp,gb,nbn,nbas,ac)
	call norma_wf(zbn,evbn,gb,nbn,nbas,ac,ispin,indg,kp)
      if(maswrk) then
        if(iprn.ge.1) then
          write(iun,'(a21,i4,3(1x,f8.4))')'Band energies at k = ',k,
     &                                    (kp(i),i=1,3)
          if(iheg) then
            do ie=1,nbas
	        iv=indg(ie)
	        v=kp+gbs(:,iv)
	        vabs=dot_product(v,v)*pi2a**2
	        write(iun,*)'e_emp e_found ',ie,vabs,ebn(ie)
            enddo   !!! over ie
          else
            do ie=1,nbn
	        write(iun,'(a17,i4,f18.13,3(1x,e12.5))')
     &	        'e d_val d_slp NRM',ie,ebn(ie),ac(1,ie),ac(2,ie),
     &                              ac(3,ie)
            enddo   !!! over ie
          endif
        else if(iprn==0.and.k==1) then
          write(iun,*)'Band energies at k = ',1
          do ie=1,nbn
	      write(iun,'(a17,i4,f18.13,3(1x,e12.5))')
     &	      'e d_val d_slp NRM',ie,ebn(ie),ac(1,ie),ac(2,ie),ac(3,ie)
          enddo   !!! over ie
	  endif
	endif
	err_val=max(err_val,maxval(ac(1,:)))
	err_slo=max(err_slo,maxval(ac(2,:)))
	deallocate(ac)
	deallocate(gb,yr,yi,ind_y)
      end
