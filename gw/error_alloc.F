      subroutine error_alloc(subr,arr,ndim,n1,n2,n3,n4,n5,n6,n7,length)
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
      integer, intent(in) :: ndim,n1,n2,n3,n4,n5,n6,n7,length
	integer :: i
	real*8 :: mem,fac
	print '(a11,i5,a11,i5,a9,i5,a10,i5,a11,i5)',' From me = ',me,
     &	' me3_tau = ',me3_tau,' me3_k = ',me3_k,' me4_kk = ',me4_kk,
     &  ' me4_pbr = ',me4_pbr
      print '(a23,a30)',' Error in subroutine : ',subr
      print '(a23,a30)',' Allocation of Array : ',arr
      fac=real(length)/1024.d0**3
      mem=1.d0
      do i=1,ndim
        if(i==1) print *,' Dimension 1 :',n1
        if(i==2) print *,' Dimension 2 :',n2
        if(i==3) print *,' Dimension 3 :',n3
        if(i==4) print *,' Dimension 4 :',n4
        if(i==5) print *,' Dimension 5 :',n5
        if(i==6) print *,' Dimension 6 :',n6
        if(i==7) print *,' Dimension 7 :',n7
        if(i==1) mem=mem*real(n1)
        if(i==2) mem=mem*real(n2)
        if(i==3) mem=mem*real(n3)
        if(i==4) mem=mem*real(n4)
        if(i==5) mem=mem*real(n5)
        if(i==6) mem=mem*real(n6)
        if(i==7) mem=mem*real(n7)
      enddo
      print '(a24,f15.6,a3)',' Total memory required :',mem*fac,' GB'
#ifdef MPI
      call mpi_abort(MPI_COMM_WORLD,ierr)
#else
      stop      
#endif
      end