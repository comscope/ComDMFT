      subroutine error_dealloc(subr,arr)
	use atom_mod
	use etot_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
#ifdef MPI
      include 'mpif.h'
	integer :: ierr
#endif
      character*30, intent(in) :: subr,arr
	print '(a11,i5,a11,i5,a9,i5,a10,i5,a11,i5)',' From me = ',me,
     &	' me3_tau = ',me3_tau,' me3_k = ',me3_k,' me4_kk = ',me4_kk,
     &  ' me4_pbr = ',me4_pbr
      print '(a23,a30)',' Error in subroutine : ',subr
      print '(a25,a30)',' Deallocation of Array : ',arr
#ifdef MPI
      call mpi_abort(MPI_COMM_WORLD,ierr)
#else
      stop      
#endif
      end