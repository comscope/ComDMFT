      subroutine bands_sig_0
c	Recalculation of the band structure with an additional 
C     self energy term
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
      implicit none
      logical :: correlation
	integer :: k,ndim,ispin,ind_tau,i,ind_k,k1,k2,n0
	complex*16, allocatable :: a(:,:,:),tmp(:,:,:),tmp1(:,:)
	correlation=.false.
      if(ubi==' qp'.or.ubi==' gw'.or.ubi=='gwg') correlation=.true.
	n0=nbndf
	allocate(tmp(n0,n0,2),a(n0,n0,npnt))
	do ispin=1,nspin
	  k1=n3_mpi_k(me3_k+1)+1
	  k2=n3_mpi_k(me3_k+1)+ndim3_k(me3_k+1)
	  do k=1,k1-1
	    e_bnd(:,k,ispin)=0.d0
	  enddo
	  do k=k2+1,npnt
	    e_bnd(:,k,ispin)=0.d0
	  enddo
	  a=(0.d0,0.d0)
	  do ind_k=1,ndim3_k(me3_k+1)
	    k=n3_mpi_k(me3_k+1)+ind_k
	    ndim=n_bnd(k,ispin)
	    tmp=(0.d0,0.d0)
	    call unpack_hermit(tmp,v_h_bb(1,1,ind_k,ispin),ndim,n0,
     &		                 n0,(0.d0,0.d0),(1.d0,0.d0))
	    do i=1,ndim
	      tmp(i,i,1)=tmp(i,i,1)+dcmplx(e_bnd(i,k,ispin),0.d0)
	    enddo
	    call unpack_hermit(tmp,v_h_bb_old(1,1,ind_k,ispin),ndim,n0,
     &		                 n0,(1.d0,0.d0),(-1.d0,0.d0))
	    call unpack_hermit(tmp,sigx_solid_k(1,1,ind_k,ispin),ndim,
     &	                     n0,n0,(1.d0,0.d0),(1.d0,0.d0))
	    if(ubi_0=='dft') then
     	      call unpack_hermit(tmp,v_xc_bb_old(1,1,ind_k,ispin),
     &                         ndim,n0,n0,(1.d0,0.d0),(-1.d0,0.d0))
     	    else
     	      call unpack_hermit(tmp,sigx_solid_k_old(1,1,ind_k,ispin),
     &                         ndim,n0,n0,(1.d0,0.d0),(-1.d0,0.d0))
          endif
c ---------------------------------------------------------------------
          a(1:ndim,1:ndim,k)=tmp(1:ndim,1:ndim,1)
c -------- Spectrum of complex matrices ---------------------------
	    if(me_tau_kk_pbr==0) call eig_val_solver(ndim,n0,a(1,1,k),
     &	                                           e_bnd(1,k,ispin))
          if(nproc_tau/=1) then
	      call brdcst(msgdbl,e_bnd(1,k,ispin),8*ndim,0,
     &	                comm_tau_kk_pbr)
	      call brdcst(msgdbl,a(1,1,k),16*n0*n0,0,comm_tau_kk_pbr)
          endif
          if(n0>ndim) e_bnd(n_bnd(k,ispin)+1:nbndf,k,ispin)=
     &                           e_bnd(n_bnd(k,ispin),k,ispin)+0.2d0
          if(iprn.ge.1.and.maswrk) then
            write(iun,'(a23,i4,3(1x,f8.4))')'X-only spectrum at k = ',k,
     &                                      (pnt(i,k),i=1,3)
            do i=1,ndim
	        write(iun,*)'evalue ',i,e_bnd(i,k,ispin)
            enddo   !!! over ie
	    endif
	  enddo  !! over ind_k    
	  if(nproc_k/=1) call dgop(a,2*n0**2*npnt,'  +',comm_pnt)
	  do ind_k=1,ndim3_k(me3_k+1)
	    k=n3_mpi_k(me3_k+1)+ind_k
	    ndim=n_bnd(k,ispin)
c ------ Recalculation the Z_BND ---------------------------------------
	    allocate(tmp1(nfun,ndim))
	    call zlacpy('a',nfun,ndim,z_bnd(1,1,ind_k,ispin),nfun,tmp1,
     &	              nfun)
	    call zgemm('n','n',nfun,ndim,ndim,(1.d0,0.d0),tmp1,nfun,
     &	             a(1,1,k),n0,(0.d0,0.d0),z_bnd(1,1,ind_k,ispin),
     &               nfun)
	    if(ubi==' qp') then
	      call zlacpy('a',nfun,ndim,z_qp(1,1,ind_k,ispin),nfun,tmp1,
     &	                nfun)
	      call zgemm('n','n',nfun,ndim,ndim,(1.d0,0.d0),tmp1,nfun,
     &	               a(1,1,k),n0,(0.d0,0.d0),z_qp(1,1,ind_k,ispin),
     &                 nfun)
	    endif
	    deallocate(tmp1)
c ------ Recalculation the EV_BND -----------------------------------
	    allocate(tmp1(nbask(k),ndim))
	    call zlacpy('a',nbask(k),ndim,ev_bnd(1,1,ind_k,ispin),nbasmpw,
     &		          tmp1,nbask(k))
	    call zgemm('n','n',nbask(k),ndim,ndim,(1.d0,0.d0),tmp1,
     &		         nbask(k),a(1,1,k),n0,(0.d0,0.d0),
     &               ev_bnd(1,1,ind_k,ispin),nbasmpw)
	    if(ubi==' qp') then
	      call zlacpy('a',nbask(k),ndim,ev_qp(1,1,ind_k,ispin),
     &	                nbasmpw,tmp1,nbask(k))
	      call zgemm('n','n',nbask(k),ndim,ndim,(1.d0,0.d0),tmp1,
     &	 	           nbask(k),a(1,1,k),n0,(0.d0,0.d0),
     &                 ev_qp(1,1,ind_k,ispin),nbasmpw)
	    endif
	    deallocate(tmp1)
	    allocate(tmp1(ndim,ndim))
c ------ Recalculation the V_HARTREE ------------------------------------
          if(ubi/='dft') then
	      call unpack_hermit(tmp,v_h_bb(1,1,ind_k,ispin),ndim,n0,
     &		                   n0,(0.d0,0.d0),(1.d0,0.d0))
	      call zgemm('n','n',ndim,ndim,ndim,(1.d0,0.d0),tmp,n0,
     &	               a(1,1,k),n0,(0.d0,0.d0),tmp1,ndim)
	      call zgemm('c','n',ndim,ndim,ndim,(1.d0,0.d0),a(1,1,k),
     &	               n0,tmp1,ndim,(0.d0,0.d0),tmp,n0)
	      call pack_hermit(tmp,v_h_bb(1,1,ind_k,ispin),ndim,n0,
     &		                 n0,0.d0,1.d0)
          endif
c ------ Recalculation the V_XC ------------------------------------
          if(ubi/='dft') then
	      call unpack_hermit(tmp,v_xc_bb(1,1,ind_k,ispin),ndim,n0,
     &		                   n0,(0.d0,0.d0),(1.d0,0.d0))
	      call zgemm('n','n',ndim,ndim,ndim,(1.d0,0.d0),tmp,n0,
     & 	               a(1,1,k),n0,(0.d0,0.d0),tmp1,ndim)
	      call zgemm('c','n',ndim,ndim,ndim,(1.d0,0.d0),a(1,1,k),
     &	               n0,tmp1,ndim,(0.d0,0.d0),tmp,n0)
	      call pack_hermit(tmp,v_xc_bb(1,1,ind_k,ispin),ndim,n0,
     &		                 n0,0.d0,1.d0)
          endif
c ------ Recalculation the SIGX_SOLID_K --------------------------------
	    call unpack_hermit(tmp,sigx_solid_k(1,1,ind_k,ispin),ndim,
     &	                     n0,n0,(0.d0,0.d0),(1.d0,0.d0))
	    call zgemm('n','n',ndim,ndim,ndim,(1.d0,0.d0),tmp,n0,
     &	             a(1,1,k),n0,(0.d0,0.d0),tmp1,ndim)
	    call zgemm('c','n',ndim,ndim,ndim,(1.d0,0.d0),a(1,1,k),
     &	             n0,tmp1,ndim,(0.d0,0.d0),tmp,n0)
	    call pack_hermit(tmp,sigx_solid_k(1,1,ind_k,ispin),ndim,
     &	                   n0,n0,0.d0,1.d0)
          if(correlation) then
            if(allocated(q_qp)) then
c ------ Recalculation the Q_QP ----------------------------------
	        call zlacpy('a',ndim,ndim,q_qp(1,1,ind_k,ispin),nbndf,
     &		              tmp1,ndim)
	        call zgemm('c','n',ndim,ndim,ndim,(1.d0,0.d0),a(1,1,k),n0,
     &	                 tmp1,ndim,(0.d0,0.d0),q_qp(1,1,ind_k,ispin),n0)
	      endif
c ------ Recalculation the SIG_C_TAU ----------------------------------
            do ind_tau=1,ndim3_tau
	        call ferm_unpack_tau(tmp,
     &			                   sig_c_tau(1,1,1,ind_tau,ind_k,ispin),
     &                             ndim,n0,n0,3)
	        do i=1,2
	          call zgemm('n','n',ndim,ndim,ndim,(1.d0,0.d0),
     &			           tmp(1,1,i),n0,a(1,1,k),n0,(0.d0,0.d0),
     &                     tmp1,ndim)
	          call zgemm('c','n',ndim,ndim,ndim,(1.d0,0.d0),a(1,1,k),
     &	                   n0,tmp1,ndim,(0.d0,0.d0),tmp(1,1,i),n0)
	        enddo
			  sig_c_tau(:,:,:,ind_tau,ind_k,ispin)=0.d0
		      call ferm_pack_tau(tmp,
     &			                 sig_c_tau(1,1,1,ind_tau,ind_k,ispin),
     &		                     ndim,n0,n0)
	      enddo
	    endif
	    deallocate(tmp1)
	  enddo   !! over ind_k  
	enddo  !! over ispin        
	if(nproc_k/=1) then
	  call dgop(e_bnd,n0*npnt*nspin,'  +',comm_pnt)
	endif
	do ispin=1,nspin
	  if(ubi==' hf') call output_eigen_states('__hf',2,ispin)
	  if(ubi==' qp') call output_eigen_states('__qp',2,ispin)
	  if(ubi==' gw') call output_eigen_states('__gw',2,ispin)
	  if(ubi=='gwg') call output_eigen_states('_gwg',2,ispin)
	enddo
	deallocate(tmp,a)
      end
