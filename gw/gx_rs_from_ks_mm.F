      subroutine gx_rs_from_ks_mm(gr,gc,ispin)
c	-----------  G_MM from G_full_0 ----------------
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer, intent(in) :: ispin
	real*8, intent(out) :: gr(maxel,maxel,ndim_nrr(me3_k+1))
	complex*16, intent(out) :: gc(maxel,maxel,ndim_nrr(me3_k+1))
	integer :: k,k0,nbnd,ind,j,ia,ib,ind_k,ig,iatom,isort,i_star,j1,
     &	         j0,ir,ir0,jatom,jsort,i0,i1,ind_ir0,i
	real*8 :: pi2,v(3),fc
	complex*16 :: cc
      complex*16, allocatable :: tmp(:,:),g_tmp(:,:),tmp2(:,:),
     &                           g_mm(:,:,:)
      allocate(g_tmp(nfun,nfun))
      allocate(tmp2(nbndf,nfun))
      allocate(tmp(nfun,nfun))
      allocate(g_mm(nqdiv,nfun,nfun))
	g_mm=(0.d0,0.d0)
	pi2=pi+pi
	do ind_k=1,ndim3_k(me3_k+1)
	  k0=n3_mpi_k(me3_k+1)+ind_k
	  nbnd=n_bnd(k0,ispin)
        g_tmp=(0.d0,0.d0)
	  if(ubi=='dft'.or.ubi==' hf') then
	    do j=1,nfun
	      do i=1,nbnd
	        tmp2(i,j)=g_full_00(i,ind_k,ispin)
     &	               *conjg(z_bnd(j,i,ind_k,ispin))
	      enddo
	    enddo
	  else
	    call zgemm('n','c',nbnd,nfun,nbnd,(1.d0,0.d0),
     &	             g_full_0(1,1,ind_k,ispin),nbndf,
     &               z_bnd(1,1,ind_k,ispin),nfun,(0.d0,0.d0),tmp2,
     &               nbndf)
        endif
	  call zgemm('n','n',nfun,nfun,nbnd,(1.d0,0.d0),
     &	           z_bnd(1,1,ind_k,ispin),nfun,tmp2,nbndf,(0.d0,0.d0),
     &             g_tmp,nfun)
	  do i_star=1,k_star(k0)
	    k=k_list(i_star,k0)
	    ig=k_group(k)
	    call zone1_number(pnt(:,k),rb0,ndiv,ind)
c ------------ Index from the left ------------------------------
	    do j=1,nfun
		    call sym_val_left(ig,g_tmp(:,j),tmp(:,j),1)
	      do iatom=1,natom
	        isort=is(iatom)
	        v=tshift(:,iatom,ig)
	        fc=pi2*dot_product(pnt(:,k),v)
	        cc=dcmplx(cos(fc),sin(fc))
	        ia=io_lem(iatom)
	        ib=ia+lfunm(isort)-1
	        tmp(ia:ib,j)=tmp(ia:ib,j)*cc
	      enddo
	    enddo
c ------------ Index from the right  ------------------------------
	    do j=1,nfun
		    call sym_val_left(ig,tmp(j,:),g_mm(ind,j,:),2)
	      do iatom=1,natom
	        isort=is(iatom)
	        v=tshift(:,iatom,ig)
	        fc=-pi2*dot_product(pnt(:,k),v)
	        cc=dcmplx(cos(fc),sin(fc))
	        ia=io_lem(iatom)
	        ib=ia+lfunm(isort)-1
	        g_mm(ind,j,ia:ib)=g_mm(ind,j,ia:ib)*cc
	      enddo
	    enddo
	  enddo
	enddo
	deallocate(tmp,g_tmp,tmp2)
	if(nproc_k/=1) call dgop(g_mm,2*nfun**2*nqdiv,'  +',comm_pnt)
c ----------------------------------------------------------
	j1=ndim_nfun(me3_k+1)
	j0=n_mpi_nfun(me3_k+1)+1
	call fft3(ndiv(1),ndiv(2),ndiv(3),nfun*j1,g_mm(1,1,j0),1)
	do ir=1,j0-1
	  g_mm(:,:,ir)=(0.d0,0.d0)
	enddo
	do ir=j0+j1,nfun
	  g_mm(:,:,ir)=(0.d0,0.d0)
	enddo
	if(nproc_k/=1) call dgop(g_mm,2*nfun**2*nqdiv,'  +',comm_pnt)
	g_mm=g_mm/dfloat(nqdiv)
c ----------------------------------------------------------
	do ind_ir0=1,ndim_nrr(me3_k+1)
	  ir0=n_mpi_nrr(me3_k+1)+ind_ir0
	  ir=nrr_red(1,ir0)
	  iatom=nrr_red(2,ir0)
	  jatom=nrr_red(3,ir0)
	  isort=is(iatom)
	  jsort=is(jatom)
	  i0=io_lem(iatom)-1
	  j0=io_lem(jatom)-1
	  if(irel<=1) then
	    do j1=1,lfunm(jsort)
	      do i1=1,lfunm(isort)
	        gr(i1,j1,ind_ir0)=g_mm(ir,i0+i1,j0+j1)
	      enddo
	    enddo
	  else if(irel==2) then
	    do j1=1,lfunm(jsort)
	      do i1=1,lfunm(isort)
	        gc(i1,j1,ind_ir0)=g_mm(ir,i0+i1,j0+j1)
	      enddo
	    enddo
	  endif
	enddo
	deallocate(g_mm)
      end
