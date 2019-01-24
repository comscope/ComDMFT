!> \brief Distribute a set of tasks over processors
!>
!> If there are N tasks and NPR processors then each processor is given a
!> number of tasks ISH=N/NPR. The first IREST processors are given one extra
!> task.
!>
      subroutine size_shift_par(n,npr,ndim_pr,n_mpi_pr)
      implicit none
	integer, intent(in) :: n   !> The number of "tasks"
	integer, intent(in) :: npr !> The number of processors
	integer, intent(out) :: ndim_pr(npr)  !> The number of tasks for each
                                              !> processor
	integer, intent(out) :: n_mpi_pr(npr) !> The first task for each
                                              !> processor
	integer :: ish,irest,ind,i
	ish=n/npr
	irest=n-ish*npr
	ndim_pr(1:irest)=ish+1
	ndim_pr(irest+1:npr)=ish
	ind=0
	do i=1,npr
	  n_mpi_pr(i)=ind
	  ind=ind+ndim_pr(i)
	enddo
	end
