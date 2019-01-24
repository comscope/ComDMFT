      subroutine v_coul_full_ii(nplwg,v_full,ind_q)
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
      implicit none
	integer, intent(in) :: nplwg,ind_q
	complex*16, intent(out) :: v_full(nplwg,nplwg)
	integer :: n_pbt,ind_j,j
      real*8, allocatable :: v0(:,:)
	n_pbt=n_pbmt+nplwg
	allocate(v0(nplwg,nplwg))
	v0=0.d0
	do ind_j=1,ndim_pbt(me3_tau+1)
	  j=n_mpi_pbt(me3_tau+1)+ind_j
	  if(j<=n_pbmt) cycle
	  if(j>n_pbt) cycle
	  v0(1:nplwg,j-n_pbmt)=v_coul_pb(n_pbmt+1:n_pbt,ind_j,ind_q)
	enddo
	if(nproc_tau/=1) call dgop(v0,nplwg**2,'  +',comm_pnk)
	call unpack_hermit(v_full,v0,nplwg,nplwg,nplwg,(0.d0,0.d0),
     &	                 (1.d0,0.d0))
      deallocate(v0)
      end