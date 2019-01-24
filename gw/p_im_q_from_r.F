      subroutine p_im_q_from_r(p_mi,ind_tau,p)
c     key = 0 : contribution to P_W_SOLID_Q_TAU -----------------
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer, intent(in) :: ind_tau
	real*8, intent(in) :: p_mi(nr_red,n_pbmt,ndim3_kk(me3_k+1))
      real*8, intent(inout) :: p(n_pbtot,n_pbtot,ndim3_tau,
     &                           maxval(ndim3_k))
	integer :: ind,i,i0,j,ia,ib,ind_k,k,ig,jatom,jsort,ind_ir,ir,k0,k1
	real*8 :: v(3),fc,pi2
	complex*16 :: cc
	complex*16, allocatable :: tmp(:,:),tmp1(:,:),tmp2(:,:),s(:,:),
     &	                         p_c(:,:),pcc(:,:,:,:)
	pi2=pi+pi
	allocate(p_c(nqdiv,n_pbmt))
	allocate(pcc(n_pbmt,nr_red,k_star_max,maxval(ndim3_k)))
	do i0=1,nr_red
	  p_c=(0.d0,0.d0)
	  do ind_ir=1,ndim3_kk(me3_k+1)
	    ir=n3_mpi_kk(me3_k+1)+ind_ir
	    do j=1,n_pbmt
	      p_c(ir,j)=p_mi(i0,j,ind_ir)
	    enddo
	  enddo
	  if(nproc_k/=1) call dgop(p_c,2*nqdiv*n_pbmt,'  +',comm_pnt)
	  call fft3(ndiv(1),ndiv(2),ndiv(3),n_pbmt,p_c,-1)
	  do ind_k=1,ndim3_k(me3_k+1)
	    k=n3_mpi_k(me3_k+1)+ind_k
     	    do k0=1,k_star(k)
	      k1=k_list(k0,k)
	      call zone1_number(pnt(1,k1),rb0,ndiv,ind)
	      pcc(:,i0,k0,ind_k)=p_c(ind,:)
	    enddo
	  enddo
	enddo  !! over i0
	deallocate(p_c)
	allocate(tmp(nr_full,n_pbmt))
	allocate(tmp1(n_pbmt,nplw_rm))
	allocate(tmp2(n_pbmt,nplw_gw))
c ---------------------------------------------------------------	
	do ind_k=1,ndim3_k(me3_k+1)
	  k=n3_mpi_k(me3_k+1)+ind_k
c ------- We fill in the first index completely -------------
	  do i=1,nr_full
	    i0=ir_ref(i)
	    ig=r_group(i)
	    ind=k_sym_0(k,ig)  ! IND=A^(-1)*K  in NQDIV-list
	    ind=ind_k_star_back(ind) ! position in a given K-star
		  call sym_vmi_left(ig,pcc(:,i0,ind,ind_k),tmp(i,:))
	    do jatom=1,natom
	      jsort=is(jatom)
	      v=tshift(:,jatom,ig)-tshift_r(:,i,ig)
	      fc=pi2*dot_product(pnt(:,k),v)
	      cc=dcmplx(cos(fc),sin(fc))
	      ia=iopb(jatom)
	      ib=ia+n_pbmt0(jsort)-1
	      tmp(i,ia:ib)=tmp(i,ia:ib)*cc
	    enddo
c ----- Multiplication with the phase factors ------------------
          v(:)=r_pnt(1,i)*r0b0(:,1)+r_pnt(2,i)*r0b0(:,2)
     &	  	                       +r_pnt(3,i)*r0b0(:,3)
	    fc=2.d0*pi*dot_product(pnt(:,k),v)
	    cc=exp(dcmplx(0.d0,fc))
	    tmp(i,1:n_pbmt)=tmp(i,1:n_pbmt)*cc
	  enddo   !! over ist
c -----------------------------------------------------------
	  call fft3(nrdiv(1),nrdiv(2),nrdiv(3),n_pbmt,tmp,1)
	  tmp=tmp/nr_full
	  tmp1=(0.d0,0.d0)
	  do i=1,nplw_r(k)
	    i0=indpw_r(i,k)
	    call zone1_number(gbs(:,i0),rbas,nrdiv,ind)
	    tmp1(1:n_pbmt,i)=tmp(ind,1:n_pbmt)
	  enddo
c ------ We form the overlap matrix -----------------------------
	  allocate(s(nplw_r(k),nplwgw(k)))
        call s_overlap_m_b(s,nplw_r(k),nplwgw(k),indpw_r(1,k),
     &                     indpw_gw(1,k))
	  call zgemm('n','n',n_pbmt,nplwgw(k),nplw_r(k),(1.d0,0.d0),
     &	           tmp1,n_pbmt,s,nplw_r(k),(0.d0,0.d0),tmp2,n_pbmt)
        tmp1(:,1:nplwgw(k))=tmp2(:,1:nplwgw(k))
        deallocate(s)
	  do j=1,nplwgw(k)
	    do i=1,n_pbmt
	      p(i,n_pbmt+j,ind_tau,ind_k)=p(i,n_pbmt+j,ind_tau,ind_k)
     &	                               +dreal(tmp1(i,j))
	      p(n_pbmt+j,i,ind_tau,ind_k)=p(n_pbmt+j,i,ind_tau,ind_k)
     &  	                           +dimag(tmp1(i,j))
	    enddo
	  enddo
c --------------------------------------------------------------
	enddo  !! over ind_k
	deallocate(tmp2)
	deallocate(tmp,tmp1,pcc)
      end
