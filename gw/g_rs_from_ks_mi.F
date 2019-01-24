      subroutine g_rs_from_ks_mi(gr,gc,it,ispin,ind_tau)
c	----------  G_IM from G_full -----------------------
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer, intent(in) :: it,ispin,ind_tau
	real*8, intent(out) :: gr(nr_red,nfun,ndim3_kk(me3_k+1),2)
	complex*16, intent(out) :: gc(nr_red,nfun,ndim3_kk(me3_k+1),2,
     &	                            nrel)
	integer :: k,k0,nbnd,ind,i,i0,ind_k,ig,j,j1,j0,iatom,isort,ia,ib,
     &	         i_tau,i_star,nbas,irl,ii,ind_ir,ir
	real*8 :: v(3),fc,pi2,tu
	complex*16 :: cc
	real*8, allocatable :: gx_tau(:)
      complex*16, allocatable :: g_tmp(:,:,:),tmp(:,:),tmp1(:,:),
     &	                         tmp0(:,:),tmp2(:,:),tmp3(:),
     &                           g_mi(:,:,:,:)
	allocate(g_mi(nqdiv,nr_red,nfun,nrel))
      allocate(tmp(nr_full,nfun),tmp1(nfun,nbasmpw))
      allocate(g_tmp(nr_full,nfun,nrel),tmp3(nfun))
	allocate(tmp0(nfun,nbndf),tmp2(nbndf,nbndf))
	pi2=pi+pi
	g_mi=(0.d0,0.d0)
	i_tau=me3_tau*ndim3_tau+ind_tau-1
	if(it==2) i_tau=n_tau-i_tau
	allocate(gx_tau(nbndf))
	do ind_k=1,ndim3_k(me3_k+1)
	  k0=n3_mpi_k(me3_k+1)+ind_k
	  nbnd=n_bnd(k0,ispin)
	  nbas=nbask(k0)/nrel
        g_tmp=(0.d0,0.d0)
        if(ubi=='dft'.or.ubi==' hf') then
	    tu=tau_mesh(i_tau)
          call g_x_tau(ispin,k0,gx_tau,tu,nbnd,chem_pot)
          do j=1,nfun
            do ib=1,nbnd
              tmp0(j,ib)=z_bnd(j,ib,ind_k,ispin)*gx_tau(ib)
            enddo
          enddo
        else
		  call ferm_unpack_tau(tmp2,g_full(1,1,1,ind_tau,ind_k,ispin),
     &	                       nbnd,nbndf,nbndf,it)
	    call zgemm('n','n',nfun,nbnd,nbnd,(1.d0,0.d0),
     &	             z_bnd(1,1,ind_k,ispin),nfun,tmp2,nbndf,(0.d0,0.d0),
     &               tmp0,nfun)
        endif
	  call zgemm('n','c',nfun,nbask(k0),nbnd,(1.d0,0.d0),tmp0,nfun,
     &	           ev_bnd(1,1,ind_k,ispin),nbasmpw,(0.d0,0.d0),tmp1,
     &             nfun)
        tmp1=tmp1/sqrt(amega)
        do irl=1,nrel
	    tmp=(0.d0,0.d0)
	    do i=1,nbas
	      i0=indgb(i,k0)
	      call zone1_number(gbs(:,i0),rbas,nrdiv,ind)
            ii=(irl-1)*nbas+i
            tmp(ind,1:nfun)=tmp1(1:nfun,ii)
	    enddo
	    call fft3(nrdiv(1),nrdiv(2),nrdiv(3),nfun,tmp,-1)
c ----- Multiplication with the phase factor ------------------
          do i=1,nr_full
            v(:)=r_pnt(1,i)*r0b0(:,1)+r_pnt(2,i)*r0b0(:,2)
     &		  	                     +r_pnt(3,i)*r0b0(:,3)
	      fc=-pi2*dot_product(pnt(:,k0),v)
	      cc=exp(dcmplx(0.d0,fc))
	      g_tmp(i,1:nfun,irl)=tmp(i,1:nfun)*cc
	    enddo
	  enddo   !! over irl
c ----------------------------------------------------------------
	  do i_star=1,k_star(k0)
	    k=k_list(i_star,k0)
	    ig=k_group(k)
	    call zone1_number(pnt(:,k),rb0,ndiv,ind)
	    do j=1,nr_red
	      j1=r_pnt_ind(j)
	      j0=ip_r(j1,ig)
            do irl=1,nrel
              if(irel<=1) then
                tmp3=g_tmp(j0,:,1)
              else if(irl==1) then
                tmp3=uj(1,ig)*g_tmp(j0,:,1)+uj(3,ig)*g_tmp(j0,:,2)
              else if(irl==2) then
                tmp3=uj(2,ig)*g_tmp(j0,:,1)+uj(4,ig)*g_tmp(j0,:,2)
              endif
		      call sym_val_left(ig,tmp3,g_mi(ind,j,:,irl),1)
	        do iatom=1,natom
	          isort=is(iatom)
	          v=tshift(:,iatom,ig)-tshift_r(:,j1,ig)
	          fc=pi2*dot_product(pnt(:,k),v)
	          cc=dcmplx(cos(fc),sin(fc))
	          ia=io_lem(iatom)
	          ib=ia+lfunm(isort)-1
	          g_mi(ind,j,ia:ib,irl)=g_mi(ind,j,ia:ib,irl)*cc
	        enddo
	      enddo
	    enddo
	  enddo
	enddo
	deallocate(g_tmp)
	deallocate(tmp0,tmp,tmp1,tmp2,tmp3)
	if(nproc_k/=1) call dgop(g_mi,2*nqdiv*nr_red*nfun*nrel,'  +',
     &	                       comm_pnt)
	call fft3(ndiv(1),ndiv(2),ndiv(3),nfun*nr_red*nrel,g_mi,1)
	g_mi=g_mi/dfloat(nqdiv)
	do ind_ir=1,ndim3_kk(me3_k+1)
	  ir=n3_mpi_kk(me3_k+1)+ind_ir
	  if(irel<=1) then
	    gr(:,:,ind_ir,it)=g_mi(ir,:,:,1)
	  else if(irel==2) then
	    do irl=1,2
	      gc(:,:,ind_ir,it,irl)=g_mi(ir,:,:,irl)
	    enddo
	  endif
	enddo
	deallocate(g_mi,gx_tau)
      end
