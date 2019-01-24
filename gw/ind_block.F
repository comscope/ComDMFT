      integer function ind_block(n_mpi,ndim_mpi,i_quest)
	use parallel_mod
      implicit none
      integer, intent(in) :: ndim_mpi(nproc),n_mpi(nproc),i_quest
	integer :: i
      ind_block=0
	do i=1,nproc
	  if(i_quest>n_mpi(i).and.i_quest<=n_mpi(i)+ndim_mpi(i)) then
	    ind_block=i
	    exit
	  endif
	enddo
      end
