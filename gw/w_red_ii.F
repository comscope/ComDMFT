      subroutine w_red_ii(wii,w_v,ind_k,npl)
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
      implicit none
	integer, intent(in) :: ind_k,npl
	real*8, intent(in) :: w_v(n_pbtot,n_pbtot)
	complex*16, intent(out) :: wii(npl,npl)
	integer :: k,n0
	complex*16, allocatable :: tmp(:,:)
	k=n3_mpi_k(me3_k+1)+ind_k
      n0=n_pbmt+npl
	allocate(tmp(n0,n0))
	call boson_unpack_tau(tmp,n0,w_v,n_pbtot,n0)
	wii=tmp(n_pbmt+1:n_pbmt+npl,n_pbmt+1:n_pbmt+npl)
	deallocate(tmp)
      end
