C*************************************************************

c      SUBROUTINE mpi_bcast_int(X,LENX,rootprocess,comm)
c      use comwann_mod
c      implicit none
c      INCLUDE 'mpif.h'
c      integer,intent(in) :: lenx,comm,rootprocess
c      integer, intent(inout) :: x(lenx)
c      integer :: mxbuff
c      integer :: i,npass,length,locx,ierror
c
c      mxbuff=8192
c      NPASS = (LENX-1)/MXBUFF + 1
c      LENGTH = MXBUFF
c      LOCX = 1
c      DO I=1,NPASS
c        IF(I.EQ.NPASS) LENGTH = LENX - MXBUFF*(NPASS-1)
c        CALL MPI_bcast(X(LOCX),LENGTH,MPI_INTEGER,
c     $    rootprocess,comm,IERROR)
c        if (ierror .ne. MPI_SUCCESS) then
c          write(179,*) 'mpi_bcast_int fail',ierror
c          call flush(179)
c          call ending
c        endif
c        LOCX = LOCX + MXBUFF
c      enddo
c      END


c      SUBROUTINE mpi_bcast_dble(X,LENX,rootprocess,comm)
c      use comwann_mod
c      implicit none
c      INCLUDE 'mpif.h'
c      integer,intent(in) :: lenx,comm,rootprocess
c      double precision, intent(inout) :: x(lenx)
c      integer :: mxbuff
c      integer :: i,npass,length,locx,ierror
c
c      mxbuff=4096
c      NPASS = (LENX-1)/MXBUFF + 1
c      LENGTH = MXBUFF
c      LOCX = 1
c      DO I=1,NPASS
c        IF(I.EQ.NPASS) LENGTH = LENX - MXBUFF*(NPASS-1)
c        CALL MPI_bcast(X(LOCX),LENGTH,MPI_double_precision,
c     $    rootprocess,comm,IERROR)
c        if (ierror .ne. MPI_SUCCESS) then
c          write(179,*) 'mpi_bcast_int fail',ierror
c          call flush(179)
c          call ending
c        endif
c        LOCX = LOCX + MXBUFF
c      enddo
c      END

c      SUBROUTINE mpi_bcast_dcmplx(X,LENX,rootprocess,comm)
c      use comwann_mod
c      implicit none
c      INCLUDE 'mpif.h'
c      integer,intent(in) :: lenx,comm,rootprocess
c      complex*16, intent(inout) :: x(lenx)
c      integer :: mxbuff
c      integer :: i,npass,length,locx,ierror
c
c      mxbuff=2048
c      NPASS = (LENX-1)/MXBUFF + 1
c      LENGTH = MXBUFF
c      LOCX = 1
c      DO I=1,NPASS
c        IF(I.EQ.NPASS) LENGTH = LENX - MXBUFF*(NPASS-1)
c        CALL MPI_bcast(X(LOCX),LENGTH,MPI_double_complex,
c     $    rootprocess,comm,IERROR)
c        if (ierror .ne. MPI_SUCCESS) then
c          write(179,*) 'mpi_bcast_int fail',ierror
c          call flush(179)
c          call ending
c        endif
c        LOCX = LOCX + MXBUFF
c      enddo
c      END

c      SUBROUTINE mpi_allreduce_int(X,LENX,OP,comm)
c      use comwann_mod
c      implicit none
c      INCLUDE 'mpif.h'
c      integer,intent(in) :: lenx,op,comm
c      integer, intent(inout) :: x(lenx)
c      integer :: mxbuff
c      integer,allocatable :: buff(:)
c      integer :: i,j,npass,length,locx,ierror
c
c      mxbuff=8192
c      allocate(buff(mxbuff))
c      NPASS = (LENX-1)/MXBUFF + 1
c      LENGTH = MXBUFF
c      LOCX = 1
c      DO I=1,NPASS
c        IF(I.EQ.NPASS) LENGTH = LENX - MXBUFF*(NPASS-1)
c        CALL MPI_ALLREDUCE(X(LOCX),BUFF,LENGTH,MPI_INTEGER,
c     $    OP,comm,IERROR)
c        if (ierror .ne. MPI_SUCCESS) then
c          write(179,*) 'mpi_allreduce_int fail',ierror
c          call flush(179)
c          call ending
c        endif
c        do j=1,length
c          x(locx+j-1)=buff(j)
c        enddo
c        LOCX = LOCX + MXBUFF
c      enddo
c      deallocate(buff)
c      END

      SUBROUTINE mpi_allreduce_dble(X,LENX,OP,comm)
      use comwann_mod
      implicit none
      INCLUDE 'mpif.h'
      integer,intent(in) :: lenx,op,comm
      double precision, intent(inout) :: x(lenx)
      integer :: mxbuff
      double precision, allocatable :: buff(:)
      integer :: i,npass,length,locx,ierror

      mxbuff=4096
      allocate(buff(mxbuff))
      NPASS = (LENX-1)/MXBUFF + 1
      LENGTH = MXBUFF
      LOCX = 1
      DO I=1,NPASS
        IF(I.EQ.NPASS) LENGTH = LENX - MXBUFF*(NPASS-1)
        CALL MPI_ALLREDUCE(X(LOCX),BUFF,LENGTH,MPI_double_precision,
     $    OP,comm,IERROR)
        if (ierror .ne. MPI_SUCCESS) then
c          write(179,*) 'mpi_allreduce_dble fail',ierror
c          call flush(179)
          call ending
        endif
        CALL DCOPY(LENGTH,BUFF,1,X(LOCX),1)
        LOCX = LOCX + MXBUFF
      enddo
      deallocate(buff)
      END

      SUBROUTINE mpi_allreduce_dcmplx(X,LENX,OP,comm)
      use comwann_mod
      implicit none
      INCLUDE 'mpif.h'
      integer,intent(in) :: lenx,op,comm
      complex*16, intent(inout) :: x(lenx)
      integer :: mxbuff
      complex*16,allocatable :: buff(:)
      integer :: i,npass,length,locx,ierror


      mxbuff=2048
      allocate(buff(mxbuff))
      NPASS = (LENX-1)/MXBUFF + 1
      LENGTH = MXBUFF
      LOCX = 1
      DO I=1,NPASS
        IF(I.EQ.NPASS) LENGTH = LENX - MXBUFF*(NPASS-1)
        CALL MPI_ALLREDUCE(X(LOCX),BUFF,LENGTH,MPI_double_complex,
     $    OP,comm,IERROR)
        if (ierror .ne. MPI_SUCCESS) then
c          write(179,*) 'mpi_allreduce_dble fail',ierror
c          call flush(179)
          call ending
        endif
        CALL ZCOPY(LENGTH,BUFF,1,X(LOCX),1)
        LOCX = LOCX + MXBUFF
      enddo
      deallocate(buff)
      END

