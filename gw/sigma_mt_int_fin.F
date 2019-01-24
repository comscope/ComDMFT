      subroutine sigma_mt_int_fin(ind_k,k,ispin)
	use atom_mod
	use manager_mod
	use solid_mod
      implicit none
      integer, intent(in) :: ind_k,k,ispin
      integer :: i,j
	complex*16, allocatable :: sig_x1(:,:),tmp(:,:)
	allocate(tmp(nfun,nbndf))
	allocate(sig_x1(nbndf,nbndf))
      call zgemm('n','n',nfun,n_bnd(k,ispin),nbask(k),(1.d0,0.d0),
     &           pvg(1,1,ind_k),nfun,ev_bnd(1,1,ind_k,ispin),nbasmpw,
     &           (0.d0,0.d0),tmp,nfun)
      call zgemm('c','n',n_bnd(k,ispin),n_bnd(k,ispin),nfun,(1.d0,0.d0),
     &           z_bnd(1,1,ind_k,ispin),nfun,tmp,nfun,(0.d0,0.d0),
     &           sig_x1,nbndf)
	deallocate(tmp)
     	do j=1,n_bnd(k,ispin)
	  do i=1,j
	    sig_x1(i,j)=sig_x1(i,j)+conjg(sig_x1(j,i))
	  enddo
	enddo
	call pack_hermit(sig_x1,sigx_solid_k(1,1,ind_k,ispin),
     &	               n_bnd(k,ispin),nbndf,nbndf,1.d0,1.d0)
	deallocate(sig_x1)
      end