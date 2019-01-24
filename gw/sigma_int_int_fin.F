      subroutine sigma_int_int_fin(ind_k,k,ispin)
	use atom_mod
	use manager_mod
	use solid_mod
      implicit none
      integer, intent(in) :: ind_k,k,ispin
	complex*16, allocatable :: tmp(:,:),tmp1(:,:)
	allocate(tmp(nbasmpw,nbndf))
	allocate(tmp1(nbndf,nbndf))
      call zgemm('n','n',nbask(k),n_bnd(k,ispin),nbask(k),(1.d0,0.d0),
     &           pgg(1,1,ind_k),nbasmpw,ev_bnd(1,1,ind_k,ispin),nbasmpw,
     &           (0.d0,0.d0),tmp,nbasmpw)
      call zgemm('c','n',n_bnd(k,ispin),n_bnd(k,ispin),nbask(k),
     &           (1.d0,0.d0),ev_bnd(1,1,ind_k,ispin),nbasmpw,tmp,
     &           nbasmpw,(0.d0,0.d0),tmp1,nbndf)
	call pack_hermit(tmp1,sigx_solid_k(1,1,ind_k,ispin),
     &	               n_bnd(k,ispin),nbndf,nbndf,1.d0,1.d0)
	deallocate(tmp,tmp1)
      end