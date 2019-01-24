      subroutine g_mt_from_g_full(correlation)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use lda_dmft_mod
      implicit none
      logical, intent(in) :: correlation
      integer :: k,iatom,isort,ie,ind_k,n,ispin,i,j
      complex*16 :: wei
      complex*16, allocatable :: tmp(:,:),tmp1(:,:),tmp2(:,:)
      g_loc_0=(0.d0,0.d0)
      allocate(tmp(maxel,nbndf))
      allocate(tmp1(nbndf,nbndf))
      allocate(tmp2(nbndf,nbndf))
      do ispin=1,nspin
        do ind_k=1,ndim3_k(me3_k+1)
          k=n3_mpi_k(me3_k+1)+ind_k
          n=n_bnd(k,ispin)
          wei=dcmplx(wgt(k),0.d0)
          tmp1=(0.d0,0.d0)
          if(.not.correlation) then
            if (is_lda_dmft) then
              do i=1,n
                do j=1,n
                  tmp1(i,j)=g_full_0(i,j,ind_k,ispin)
                enddo
              enddo
            else
              do i=1,n
                tmp1(i,i)=g_full_00(i,ind_k,ispin)
              enddo
            endif
          else
            tmp1(1:n,1:n)=g_full_0(1:n,1:n,ind_k,ispin)
          endif
          do iatom=1,natom
            isort=is(iatom)
            ie=io_lem(iatom)
            call zgemm('n','n',lfunm(isort),n,n,(1.d0,0.d0),
     &        z_bnd(ie,1,ind_k,ispin),nfun,tmp1,nbndf,
     &        (0.d0,0.d0),tmp,maxel)
            call zgemm('n','c',lfunm(isort),lfunm(isort),n,wei,tmp,
     &        maxel,z_bnd(ie,1,ind_k,ispin),nfun,(1.d0,0.d0),
     &        g_loc_0(1,1,iatom,ispin),maxel)
          enddo                 !! over iatom
        enddo                   !! over ind_k
      enddo                     !! over ispin
      deallocate(tmp,tmp1,tmp2)
      if(nproc_k/=1) call DGOP(g_loc_0,2*maxel**2*natom*nspin_0,'  +',
     &  comm_pnt)
c     ------- Symmetrization -----------------------------------------------
      do ispin=1,nspin
        call symg_loc(g_loc_0(1,1,1,ispin))
        g_loc_0(:,:,:,ispin)=conjg(g_loc_0(:,:,:,ispin))
      enddo
      end
