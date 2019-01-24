      subroutine sigma_mt_mt_fin(ind_k,k,ispin)
	use atom_mod
	use manager_mod
	use solid_mod
      implicit none
      integer, intent(in) :: ind_k,k,ispin
	complex*16, allocatable :: tmp(:,:),tmp1(:,:)
	allocate(tmp(nfun,nbndf))
	allocate(tmp1(nbndf,nbndf))
	call zhemm('l','l',nfun,n_bnd(k,ispin),(1.d0,0.d0),pvv(1,1,ind_k),
     &	         nfun,z_bnd(1,1,ind_k,ispin),nfun,(0.d0,0.d0),tmp,nfun)
      call zgemm('c','n',n_bnd(k,ispin),n_bnd(k,ispin),nfun,(1.d0,0.d0),
     &           z_bnd(1,1,ind_k,ispin),nfun,tmp,nfun,(0.d0,0.d0),tmp1,
     &           nbndf)
	call pack_hermit(tmp1,sigx_solid_k(1,1,ind_k,ispin),
     &	               n_bnd(k,ispin),nbndf,nbndf,1.d0,1.d0)
	deallocate(tmp,tmp1)
      end