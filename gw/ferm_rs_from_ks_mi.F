      subroutine ferm_rs_from_ks_mi(nbm,gk,zk,ak,gr,gc,n_bn,it)
c	----------  G_IM from G_full -----------------------
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer, intent(in) :: nbm,n_bn(npnt),it
	real*8, intent(in) :: gk(nbm,nbm,maxval(ndim3_k))
      complex*16, intent(in) :: zk(nfun,nbm,npnt),ak(nbasmpw,nbm,npnt)
	real*8, intent(out) :: gr(nqdiv,nr_red,nfun,2)
	complex*16, intent(out) :: gc(nqdiv,nr_red,nfun,2,nrel)
	integer :: k,k0,nbnd,ind,i,i0,ind_k,ig,j,j1,j0,iatom,isort,ia,ib,
     &	         i_star,nbas,irl,ii
	real*8 :: v(3),fc,pi2
	complex*16 :: cc
      complex*16, allocatable :: g_tmp(:,:,:),tmp(:,:),tmp1(:,:),
     &	                         tmp0(:,:),tmp2(:,:),tmp3(:),
     &                           g_mi(:,:,:,:)
	allocate(g_mi(nqdiv,nr_red,nfun,nrel))
      allocate(tmp(nr_full,nfun),tmp1(nfun,nbasmpw))
      allocate(g_tmp(nr_full,nfun,nrel),tmp3(nfun))
	allocate(tmp0(nfun,nbm),tmp2(nbm,nbm))
	pi2=pi+pi
	g_mi=(0.d0,0.d0)
	do ind_k=1,ndim3_k(me3_k+1)
	  k0=n3_mpi_k(me3_k+1)+ind_k
	  nbnd=n_bn(k0)
	  nbas=nbask(k0)/nrel
        g_tmp=(0.d0,0.d0)
		call unpack_hermit(tmp2,gk(1,1,ind_k),nbnd,nbm,nbm,(0.d0,0.d0),
     &		               (1.d0,0.d0))
	  call zgemm('n','n',nfun,nbnd,nbnd,(1.d0,0.d0),zk(1,1,k0),nfun,
     &	           tmp2,nbm,(0.d0,0.d0),tmp0,nfun)
	  call zgemm('n','c',nfun,nbask(k0),nbnd,(1.d0,0.d0),tmp0,nfun,
     &	           ak(1,1,k0),nbasmpw,(0.d0,0.d0),tmp1,nfun)
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
	call fft3(ndiv(1),ndiv(2),ndiv(3),nfun*nr_red*nrel,g_mi,1)
	g_mi=g_mi/dfloat(nqdiv)
	if(irel<=1) then
	  gr(:,:,:,it)=g_mi(:,:,:,1)
	else if(irel==2) then
	  do irl=1,2
	    gc(:,:,:,it,irl)=g_mi(:,:,:,irl)
	  enddo
	endif
	deallocate(g_mi)
      end
