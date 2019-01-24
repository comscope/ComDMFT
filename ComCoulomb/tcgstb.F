C*************************************************************

      SUBROUTINE mpi_allreduce_int(X,LENX,OP,comm)
      use parallel_mod
      implicit none
      INCLUDE 'mpif.h'
      integer,intent(in) :: lenx,op,comm
      integer, intent(inout) :: x(lenx)
      integer :: ierror
      call mpi_allreduce(mpi_in_place,x,lenx,mpi_integer,op,comm,ierror)
      if (ierror .ne. MPI_SUCCESS) then
c        write(179,*) 'mpi_allreduce_int fail',ierror
c        call flush(179)
        call ending
      endif
      END

      SUBROUTINE mpi_allreduce_dble(X,LENX,OP,comm)
      use parallel_mod
      implicit none
      INCLUDE 'mpif.h'
      integer,intent(in) :: lenx,op,comm
      double precision, intent(inout) :: x(lenx)
      integer :: ierror
      call mpi_allreduce(mpi_in_place,x,lenx,mpi_double_precision,
     $  op,comm,ierror)
      if (ierror .ne. MPI_SUCCESS) then
c        write(179,*) 'mpi_allreduce_int fail',ierror
c        call flush(179)
        call ending
      endif
      END


      SUBROUTINE mpi_allreduce_dcmplx(X,LENX,OP,comm)
      use parallel_mod
      implicit none
      INCLUDE 'mpif.h'
      integer,intent(in) :: lenx,op,comm
      complex*16, intent(inout) :: x(lenx)
      integer :: ierror
      call mpi_allreduce(mpi_in_place,x,lenx,mpi_double_complex,
     $  op,comm,ierror)
      if (ierror .ne. MPI_SUCCESS) then
c        write(179,*) 'mpi_allreduce_int fail',ierror
c        call flush(179)
        call ending
      endif
      END            

