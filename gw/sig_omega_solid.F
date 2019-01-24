      subroutine sig_omega_solid
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
      implicit none
	integer :: k,ispin,n,ind_k,n0
	real*8, allocatable :: st(:,:,:,:)
	n0=nbndf
	allocate(st(nbndf,nbndf,2,ndim3_tau))
	do ispin=1,nspin
	  do ind_k=1,ndim3_k(me3_k+1)
	    k=n3_mpi_k(me3_k+1)+ind_k
	    n=n_bnd(k,ispin)
	    st=sig_c_tau(:,:,:,:,ind_k,ispin)
		  call from_tau_to_omega_baa(sig_c_omega(1,1,1,1,ind_k,ispin),
     &	                             n0,st,n0,n)
        enddo   !! over ind_k
      enddo  !! over ispin
      deallocate(st)
      end
