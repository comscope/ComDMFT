      subroutine w_red(w_full,ind_tau,ind_k)
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
      implicit none
	integer, intent(in) :: ind_k,ind_tau
	complex*16, intent(out) :: w_full(n_pbtot,n_pbtot)
	integer :: k,n
	k=n3_mpi_k(me3_k+1)+ind_k
      n=n_pbmt+nplwgw(k)
	call boson_unpack_tau(w_full,n_pbtot,
     &	                    p_w_solid_q_tau(1,1,ind_tau,ind_k),
     &	                    n_pbtot,n)
      end
