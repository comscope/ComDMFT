      subroutine v_coul_full_mi(nplwg,v_full,ind_q)
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
      implicit none
	integer, intent(in) :: nplwg,ind_q
	complex*16, intent(out) :: v_full(n_pbmt,nplwg)
	integer :: n_pbt,ind_j,j,i
      real*8, allocatable :: v0(:,:)
	n_pbt=n_pbmt+nplwg
	allocate(v0(n_pbt,n_pbt))
	v0=0.d0
	do ind_j=1,ndim_pbt(me3_tau+1)
	  j=n_mpi_pbt(me3_tau+1)+ind_j
	  if(j>n_pbt) cycle
	  v0(:,j)=v_coul_pb(1:n_pbt,ind_j,ind_q)
	enddo
	if(nproc_tau/=1) call dgop(v0,n_pbt**2,'  +',comm_pnk)
	do j=1,nplwg
	  do i=1,n_pbmt
	    v_full(i,j)=dcmplx(v0(i,n_pbmt+j),v0(n_pbmt+j,i))
	  enddo
	enddo
      deallocate(v0)
      end