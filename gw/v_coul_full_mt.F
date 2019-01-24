      subroutine v_coul_full_mt(v_full,ind_q)
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
      implicit none
	integer, intent(in) :: ind_q
	complex*16, intent(out) :: v_full(n_pbmt,n_pbmt)
	integer :: ind_j,j
      real*8, allocatable :: v0(:,:)
	allocate(v0(n_pbmt,n_pbmt))
	v0=0.d0
	do ind_j=1,ndim_pbt(me3_tau+1)
	  j=n_mpi_pbt(me3_tau+1)+ind_j
	  if(j>n_pbmt) cycle
	  v0(:,j)=v_coul_pb(1:n_pbmt,ind_j,ind_q)
	enddo
	if(nproc_tau/=1) call dgop(v0,n_pbmt**2,'  +',comm_pnk)
	call unpack_hermit(v_full,v0,n_pbmt,n_pbmt,n_pbmt,(0.d0,0.d0),
     &	                 (1.d0,0.d0))
      deallocate(v0)
      end