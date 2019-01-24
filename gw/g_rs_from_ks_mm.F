      subroutine g_rs_from_ks_mm(gr,gc,it,ind_tau,ispin)
c      -----------  G_MM from G_full in (t;b-t) form ----------------
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: it,ispin,ind_tau
      real*8, intent(out) :: gr(maxel,maxel,ndim_nrr(me3_k+1))
      complex*16, intent(out) :: gc(maxel,maxel,ndim_nrr(me3_k+1))
      integer :: k,k0,nbnd,ind,j,ia,ib,ind_k,ig,iatom,isort,i_tau,j1,j0,
     &           ir,ir0,jatom,jsort,i0,i1,ind_ir0,ind_i,i,m0,m,jj0,jj
      real*8 :: pi2,tu,tt(3),phase
      real*8, allocatable :: gx_tau(:)
      complex*16, allocatable :: tmp(:,:),g_tmp(:,:,:),tmp1(:,:),
     &                           tmp2(:,:),g_mm(:,:,:),v0(:,:),
     &                           g1(:,:,:),tr(:),sc(:,:,:)
      allocate(g_tmp(nfun,nfun,maxval(ndim3_k)))
      allocate(tmp2(nbndf,nfun),tmp1(nbndf,nbndf))
      pi2=pi+pi
      i_tau=me3_tau*ndim3_tau+ind_tau-1
      if(it==2) i_tau=n_tau-i_tau
      allocate(gx_tau(nbndf))
      do ind_k=1,ndim3_k(me3_k+1)
        k0=n3_mpi_k(me3_k+1)+ind_k
        nbnd=n_bnd(k0,ispin)
        if(ubi=='dft'.or.ubi==' hf') then
          tu=tau_mesh(i_tau)
          call g_x_tau(ispin,k0,gx_tau,tu,nbnd,chem_pot)
          do j=1,nfun
            do ib=1,nbnd
              tmp2(ib,j)=gx_tau(ib)*conjg(z_bnd(j,ib,ind_k,ispin))
            enddo
          enddo
        else    
      	  call ferm_unpack_tau(tmp1,g_full(1,1,1,ind_tau,ind_k,ispin),
     &                             nbnd,nbndf,nbndf,it)
          call zgemm('n','c',nbnd,nfun,nbnd,(1.d0,0.d0),tmp1,nbndf,
     &               z_bnd(1,1,ind_k,ispin),nfun,(0.d0,0.d0),tmp2,nbndf)
        endif
        call zgemm('n','n',nfun,nfun,nbnd,(1.d0,0.d0),
     &             z_bnd(1,1,ind_k,ispin),nfun,tmp2,nbndf,(0.d0,0.d0),
     &             g_tmp(1,1,ind_k),nfun)
      enddo   !! over ind_k
      deallocate(tmp2,tmp1,gx_tau)
      allocate(v0(nfun,npnt))
      allocate(g1(ndim_nfun(me3_k+1),nfun,npnt))
      do j=1,nfun
        v0=(0.d0,0.d0)
        do ind_k=1,ndim3_k(me3_k+1)
          k0=n3_mpi_k(me3_k+1)+ind_k
          v0(:,k0)=g_tmp(:,j,ind_k)
        enddo
        if(nproc_k/=1) call dgop(v0,2*npnt*nfun,'  +',comm_pnt)
        do ind_i=1,ndim_nfun(me3_k+1)
          i=n_mpi_nfun(me3_k+1)+ind_i
          g1(ind_i,j,:)=v0(i,:)
        enddo
      enddo
      deallocate(v0,g_tmp)
c ------------------------------------------------------------------
      allocate(tmp(nfun,nfun))
      allocate(tr(natom))
      allocate(g_mm(nqdiv,nfun,ndim_nfun(me3_k+1)))
      do k=1,nqdiv
        k0=i_kref(k)
        ig=k_group(k)
        call zone1_number(pnt(:,k),rb0,ndiv,ind)
        do iatom=1,natom
          tt=tshift(:,iatom,ig)
          phase=pi2*dot_product(pnt(:,k),tt)
          tr(iatom)=dcmplx(cos(phase),sin(phase))
        enddo
        tmp=(0.d0,0.d0)
        do ind_i=1,ndim_nfun(me3_k+1)
          i=n_mpi_nfun(me3_k+1)+ind_i
c ------------ Index from the right  ------------------------------
      	  call sym_val_left(ig,g1(ind_i,:,k0),tmp(:,i),2)
          do iatom=1,natom
            isort=is(iatom)
            ia=io_lem(iatom)
            ib=ia+lfunm(isort)-1
            tmp(ia:ib,i)=tmp(ia:ib,i)*conjg(tr(iatom))
          enddo
        enddo
        if(nproc_k/=1) call dgop(tmp,2*nfun**2,'  +',comm_pnt)
        tmp=transpose(tmp)
        do ind_i=1,ndim_nfun(me3_k+1)
          i=n_mpi_nfun(me3_k+1)+ind_i
c ------------ Index from the left ------------------------------
      	  call sym_val_left(ig,tmp(:,i),g_mm(ind,:,ind_i),1)
          do iatom=1,natom
            isort=is(iatom)
            ia=io_lem(iatom)
            ib=ia+lfunm(isort)-1
            g_mm(ind,ia:ib,ind_i)=g_mm(ind,ia:ib,ind_i)*tr(iatom)
          enddo
        enddo
      enddo
      deallocate(tmp,g1,tr)
c ----------------------------------------------------------
      do ind_i=1,ndim_nfun(me3_k+1)
        call fft3(ndiv(1),ndiv(2),ndiv(3),nfun,g_mm(1,1,ind_i),1)
      enddo
      g_mm=g_mm/dfloat(nqdiv)
c ----------------------------------------------------------
      allocate(sc(maxel,maxel,nrr_max))
      sc=(0.d0,0.d0)
      m0=n_mpi_nfun(me3_k+1)+1
      m=m0+ndim_nfun(me3_k+1)-1
      do ir0=1,nrr_max
        ir=nrr_red(1,ir0)
        iatom=nrr_red(2,ir0)
        jatom=nrr_red(3,ir0)
        jsort=is(jatom)
        j0=io_lem(jatom)
        if(j0>m) cycle
        j=j0+lfunm(jsort)-1
        if(j<m0) cycle
        j0=max(j0,m0)
        j=min(j,m)
        jj0=j0-n_mpi_nfun(me3_k+1)
        jj=j-n_mpi_nfun(me3_k+1)
        j0=j0-io_lem(jatom)+1
        j=j-io_lem(jatom)+1
        isort=is(iatom)
        i0=io_lem(iatom)
        i=i0+lfunm(isort)-1
        sc(1:lfunm(isort),j0:j,ir0)=g_mm(ir,i0:i,jj0:jj)
      enddo
      deallocate(g_mm)
      if(nproc_k/=1) call dgop(sc,2*maxel**2*nrr_max,'  +',comm_pnt)
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
              gr(i1,j1,ind_ir0)=sc(i1,j1,ir0)
            enddo
          enddo
        else if(irel==2) then
          do j1=1,lfunm(jsort)
            do i1=1,lfunm(isort)
              gc(i1,j1,ind_ir0)=sc(i1,j1,ir0)
            enddo
          enddo
        endif
      enddo
      deallocate(sc)
      end
