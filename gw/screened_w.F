      subroutine screened_w(w_lda)
      use atom_mod
      use manager_mod
	use solid_mod
	use parallel_mod
	implicit none
	logical, intent(in) :: w_lda
	integer :: ind_k,k,n_pbt
	complex*16, allocatable :: tmp(:,:)
	if(itest/=0) then
	  allocate(tmp(n_pbtot,n_pbtot))
	  if(maswrk) write(iun,*)' Symmetry check for P:'
	  do ind_k=1,ndim3_k(me3_k+1)
	    k=n3_mpi_k(me3_k+1)+ind_k
	    n_pbt=n_pbmt+nplwgw(k)
	    call boson_unpack_tau(tmp,n_pbtot,
     &                          p_w_solid_q_tau(1,1,1,ind_k),
     &                          n_pbtot,n_pbt)
	    call check_w_sym(k,tmp,n_pbtot,2)
        enddo
      endif
      call ws_k_nu_solid(w_lda)
	if(itest/=0) then
	  if(maswrk) write(iun,*)' Symmetry check for W:'
	  do ind_k=1,ndim3_k(me3_k+1)
	    k=n3_mpi_k(me3_k+1)+ind_k
	    n_pbt=n_pbmt+nplwgw(k)
	    call boson_unpack_tau(tmp,n_pbtot,
     &                          p_w_solid_q_tau(1,1,1,ind_k),
     &                          n_pbtot,n_pbt)
	    call check_w_sym(k,tmp,n_pbtot,2)
        enddo
      endif
      if(itest/=0) deallocate(tmp)
      call timel('**** SCREENED_W finished ***********')
      end
