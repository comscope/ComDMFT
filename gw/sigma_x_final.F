      subroutine sigma_x_final(ndim,sig_x0,ztmp,nbnd,ind_k,ispin)
	use atom_mod
	use manager_mod
	use solid_mod
      implicit none
      integer, intent(in) :: nbnd,ndim,ind_k,ispin
      complex*16, intent(in) :: sig_x0(ndim,ndim),ztmp(ndim,ndim)
	complex*16, allocatable :: tmp(:,:),tmp1(:,:)
	allocate(tmp(ndim,nbnd))
	allocate(tmp1(nbnd,nbnd))
	call zhemm('l','l',ndim,nbnd,(1.d0,0.d0),sig_x0,ndim,
     &	       ztmp,ndim,(0.d0,0.d0),tmp,ndim)
      call zgemm('c','n',nbnd,nbnd,ndim,(1.d0,0.d0),ztmp,ndim,tmp,
     &	       ndim,(0.d0,0.d0),tmp1,nbnd)
	call pack_hermit(tmp1,sigx_solid_k(1,1,ind_k,ispin),nbnd,nbnd,
     &	             nbndf,0.d0,1.d0)
	deallocate(tmp,tmp1)
      end
