      subroutine symmetrization_z      
	use atom_mod
	use manager_mod
	use parallel_mod
	use units_mod
	use solid_mod
      implicit none
	integer :: ig,ii,ispin,ind_k,k
	complex*16, allocatable :: tmp(:,:),tmp1(:,:)
	allocate(tmp(nfun,nbndf))
	allocate(tmp1(nfun,nbndf))
	do ispin=1,nspin
	  do ind_k=1,ndim3_k(me3_k+1)
	    k=n3_mpi_k(me3_k+1)+ind_k
	    tmp1=(0.d0,0.d0)
	    do ii=1,num_sym_k(k)
	      ig=list_sym_k(ii,k)
	      call sym_z_0(tmp,k,z_bnd(1,1,ind_k,ispin),n_bnd(k,ispin),ig,
     &	                 pnt(1,k))
            tmp1=tmp1+tmp
          enddo
          z_bnd(:,:,ind_k,ispin)=tmp1/num_sym_k(k)
        enddo
      enddo
	deallocate(tmp,tmp1)
      end
