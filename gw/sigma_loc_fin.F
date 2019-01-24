      subroutine sigma_loc_fin(ind_k,k,ispin)
	use atom_mod
	use manager_mod
	use solid_mod
      implicit none
      integer, intent(in) :: ind_k,k,ispin
      integer :: iatom,isort,ind,n,n0
	complex*16, allocatable :: tmp(:,:),tmp1(:,:)
	n0=nbndf
	n=n_bnd(k,ispin)
	allocate(tmp(nfun,n0))
	allocate(tmp1(n0,n0))
c	pv_x=conjg(pv_x)
      do iatom=1,natom
        isort=is(iatom)
        ind=io_lem(iatom)
	  call zhemm('l','l',lfunm(isort),n,(1.d0,0.d0),pv_x(1,1,iatom),
     &	           maxel,z_bnd(ind,1,ind_k,ispin),nfun,(0.d0,0.d0),
     &             tmp(ind,1),nfun)
      enddo   !!! over iatom
      call zgemm('c','n',n,n,nfun,(1.d0,0.d0),z_bnd(1,1,ind_k,ispin),
     &           nfun,tmp,nfun,(0.d0,0.d0),tmp1,n0)
	call pack_hermit(tmp1,sigx_solid_k(1,1,ind_k,ispin),n,n0,n0,1.d0,
     &	               1.d0)
	deallocate(tmp,tmp1)
      end