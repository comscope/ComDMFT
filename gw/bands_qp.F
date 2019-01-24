      subroutine bands_qp(x)
c     x - current chemical potential
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
      implicit none
      real*8, intent(in) :: x
      integer :: ispin,k,ndim,j,ind_k,i,il
	complex*16, allocatable :: a(:,:),b(:,:),cc(:,:)
c ------------------------------------------------------------
	allocate(a(nbndf,nbndf),b(nbndf,nbndf))
      allocate(cc(nbndf,nbndf))
      e_qp=0.d0
      ev_qp=(0.d0,0.d0)
      z_qp=(0.d0,0.d0)
      qp_spectr_weight=0.d0
	do ispin=1,nspin
	  do ind_k=1,ndim3_k(me3_k+1)
	    k=n3_mpi_k(me3_k+1)+ind_k
	    ndim=n_bnd(k,ispin)
	    call unpack_hermit(a,sig_c_0(1,1,ind_k,ispin),ndim,nbndf,
     &	                     nbndf,(0.d0,0.d0),(1.d0,0.d0))
	    call unpack_hermit(b,z_ren(1,1,ind_k,ispin),ndim,nbndf,
     &	                     nbndf,(0.d0,0.d0),(1.d0,0.d0))
c  Matrix b is inverse Z factor now - we diagonalize it to get Z^(1/2)
c ---------- Eigen pairs of Z^(-1) ----------------------------------
	    do j=1,ndim
	      do i=1,ndim
	        cc(i,j)=b(i,j)
	      enddo
	    enddo
          if(me_tau_kk_pbr==0) call eig_val_solver(ndim,nbndf,cc,
     &                                             e_qp(1,k,ispin))
          if(nproc_tau/=1) then
	      call brdcst(msgdbl,e_qp(1,k,ispin),8*ndim,0,
     &	                comm_tau_kk_pbr)
	      call brdcst(msgdbl,cc,16*nbndf**2,0,comm_tau_kk_pbr)
          endif
          if(e_qp(1,k,ispin)<=0.d0) then
            print *,'me : ',me,' problem in BANDS_QP'
            print *,' k diag(1)',k,e_qp(1,k,ispin)
          endif
          if(nbndf>ndim) e_qp(ndim+1:nbndf,k,ispin)=
     &                     e_qp(ndim,k,ispin)+0.2d0
	    b=(0.d0,0.d0)
	    do j=1,ndim
	      do i=1,ndim
	        do il=1,ndim
	          b(i,j)=b(i,j)+cc(i,il)/sqrt(e_qp(il,k,ispin))*
     &                          conjg(cc(j,il))
	        enddo
	      enddo
	    enddo
c ----- Now B = Z^(1/2) ------------------------------------------     
	    do j=1,ndim
	      a(j,j)=a(j,j)+e_bnd(j,k,ispin)-x
	    enddo
 	    call zgemm('n','n',ndim,ndim,ndim,(1.d0,0.d0),a,nbndf,b,
     &	             nbndf,(0.d0,0.d0),cc,nbndf)
	    call zgemm('n','n',ndim,ndim,ndim,(1.d0,0.d0),b,nbndf,cc,
     &	             nbndf,(0.d0,0.d0),a,nbndf)
	    do j=1,ndim
	      do i=1,ndim
	        cc(i,j)=a(i,j)
	      enddo
	    enddo
          if(me_tau_kk_pbr==0) call eig_val_solver(ndim,nbndf,cc,
     &                                             e_qp(1,k,ispin))
          if(nproc_tau/=1) then
	      call brdcst(msgdbl,e_qp(1,k,ispin),8*ndim,0,
     &	                comm_tau_kk_pbr)
	      call brdcst(msgdbl,cc,16*nbndf**2,0,comm_tau_kk_pbr)
          endif
	    do j=1,ndim
	      e_qp(j,k,ispin)=e_qp(j,k,ispin)+x
	    enddo
c --- We put renormalization Z = 1 here ------------------------	    
	    q_qp(1:ndim,1:ndim,ind_k,ispin)=cc(1:ndim,1:ndim)
c ---------------------------------------------------------------	    
	    call zgemm('n','n',nfun,ndim,ndim,(1.d0,0.d0),
     &		         z_bnd(1,1,ind_k,ispin),nfun,q_qp(1,1,ind_k,ispin),
     &               nbndf,(0.d0,0.d0),z_qp(1,1,ind_k,ispin),nfun)
	    call zgemm('n','n',nbask(k),ndim,ndim,(1.d0,0.d0),
     &		         ev_bnd(1,1,ind_k,ispin),nbasmpw,
     &               q_qp(1,1,ind_k,ispin),nbndf,(0.d0,0.d0),
     &               ev_qp(1,1,ind_k,ispin),nbasmpw)
          do i=1,ndim
            qp_spectr_weight(i,ind_k,ispin)=
     &          sum(abs(q_qp(1:ndim,i,ind_k,ispin))**2)
          enddo
	  enddo   ! over ind_k
	enddo  !! over ispin
	if(nproc_k/=1) call dgop(e_qp,nbndf*npnt*nspin,'  +',comm_pnt)
	deallocate(a,b,cc)
	do ispin=1,nspin
	  if(ubi==' qp') call output_eigen_states('__qp',3,ispin)
	  if(ubi==' gw') call output_eigen_states('__gw',3,ispin)
	  if(ubi=='gwg') call output_eigen_states('_gwg',3,ispin)
	enddo
      end
