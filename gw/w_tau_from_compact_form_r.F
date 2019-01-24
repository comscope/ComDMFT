      subroutine w_tau_from_compact_form_r(w_store,wtau,n0,n1,n)
      use atom_mod
      use manager_mod
      use parallel_mod
	use solid_mod
      implicit none
      integer, intent(in) :: n0,n,n1
      real*8, intent(in) :: w_store(n0,n0,ndim3_tau)
      real*8, intent(out) :: wtau(n1,n1,0:n_tau)
      integer :: ind_tau,i,j,i_tau
      wtau=0.d0
      do ind_tau=1,ndim3_tau
	  i_tau=me3_tau*ndim3_tau+ind_tau-1
	  do j=1,n
	    do i=1,n
            wtau(i,j,i_tau)=w_store(i,j,ind_tau)
            wtau(i,j,n_tau-i_tau)=w_store(i,j,ind_tau)
          enddo
        enddo
      enddo
	if(nproc_tau/=1) call dgop(wtau,n1**2*(n_tau+1),'  +',comm_pnk)
      end