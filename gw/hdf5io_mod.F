c-----------------------------------------------------------------------
c
      module hdf5io_mod
#ifdef USE_HDF5
        use hdf5
        use iso_c_binding
        implicit none
        private
        public :: finit, fopen, fclose, ffinal, fwrite, fread, fexist
        include "mpif.h"
c
c Module providing facilities for I/O through HDF5
c
c HDF5 provides a convenient interface to read and write data from and
c to files. In addition the file format is self documenting so that
c the files can be used with third party tools. This module requires
c an HDF5 library configured with the following options:
c
c     --enable-fortran --enable-fortran2003 --enable-parallel
c
c Note that the program uses different data distributions for different
c data sets. For the most part the code relies on replicated data
c parallelization, some data partially distributed data (the data is
c distributed over sub-sets of cores but replicated within each
c sub-set), and some data is fully distributed. Examples of these
c variables are given below:
c
c  --------------------------------------
c  replicated | partially   | distributed
c             | distributed | 
c  --------------------------------------
c  g_loc_0    | g_full_0    | g_full
c  ro_core    | g_full_00   |
c             | z_ren       |
c
c To deal with the I/O in a parallel setting this module takes the
c following approach:
c - fopen, fclose - deal with HDF5 files using MPI_COMM_WORLD
c - fread, fwrite - deal with reading and writing data assuming that the
c                   data is replicated within a given communicator
c The behavior of the read and write routines is based on the assumption
c that these routines are called collectively across the world
c communicator. The communicator passed to these routines is the one
c within which the data is replicated. At the moment rank zero of
c each communicator writes the data to file. All ranks issue the read
c call when reading the data from file to ensure the data is replicated
c in the right way upon return (we may replace this with a read and
c broadcast approach if needed). For data that is fully distributed
c "dummy" communicators of size 1 have to be provided (we may change
c that later).
c
c For example, consider a calculation that is running on 4 processors,
c and that uses 3 arrays:
c - array4(4) - is replicated on all 4 processors
c - array2(2) - is distributed over 2 sets of 2 processors each
c - array1(1) - is full distributed 
c Assume that the aggregates of array4, array2, and array1 across all
c processors happen to hold the same data. 
c Likewise we have 3 different levels of communicators
c - comm4 == MPI_COMM_WORLD
c - comm2 - communicator with 2 processors
c - comm1 - communicator with 1 processors
c Below is what the communicators and memory contents on the 4 
c processors looks like:
c
c            -------------------------------------------------
c comm4 rank |     0     |     1     |     2     |     3     |
c comm2 rank |     0     |     1     |     0     |     1     |
c comm1 rank |     0     |     0     |     0     |     0     |
c            |-----------|-----------|-----------|-----------|
c array4     | /1,2,3,4/ | /1,2,3,4/ | /1,2,3,4/ | /1,2,3,4/ |
c array2     |   /1,2/   |   /1,2/   |   /3,4/   |   /3,4/   |
c array1     |    /1/    |    /2/    |    /3/    |    /4/    |
c            -------------------------------------------------
c
c Now we have the following scenarios for the fwrite routines:
c - comm4,array4:
c   - rank 0 of comm4 writes array4 to disk
c   - this means rank 0 writes /1,2,3,4/, ranks 1 to 3 do nothing
c - comm2,array2:
c   - every rank 0 of comm2 writes its array2 to disk
c   - this means comm4 rank 0 writes /1,2/, rank 2 writes /3,4/, ranks
c     1 and 3 do nothing
c - comm1,array1:
c   - every rank saves its entire copy of array1.
c For the fread routines we have:
c - comm4,array4:
c   - every rank reads all of array4
c   - so rank 0 reads /1,2,3,4/, rank 1 read /1,2,3,4/, etc.
c - comm2,array2:
c   - every sub-set works out which slice of array2 they need
c   - all ranks in comm2 read the whole slice they are responsible for
c   - so comm4 rank 0 reads /1,2/, rank 1 reads /1,2/, rank 2 reads
c     /3,4/, and rank 3 reads /3,4/
c - comm1,array1:
c   - every sub-set works out which slice of array1 they need
c   - each rank reads the slice of array1 it needs
c   - so comm4 rank 0 reads /1/, rank 1 reads /2/, etc.
c
c In order to execute the above approach the reading and writing 
c routines need to be given a number of data items. In general will
c assume that the aggregated data structure will be an N-dimensional
c rectangular array. Each processor will hold an N-dimensional block
c of this array. To locate the position of this local block in the 
c global array an offset vector needs to be given. We also need the
c communicator across which the local block is replicated.
c So if we have a global array with dimensions NG such that 
c
c   global(1:ng(1),1:ng(2),...,1:ng(nd))
c
c a given processor will hold a local part
c
c   local(1:nl(1),1:nl(2),...,1:nl(nd))
c
c such that
c
c    global(nv(1)+1:nv(1)+nl(1),...,nv(nd)+1:nv(nd)+nl(nd))
c    == local(1:nl(1),...,1:nl(nd))
c
c Hence the read and write routines need to be given NG, NV, NL, ND
c as well as LOCAL. The GLOBAL array is the entity that lives on the
c HDF5 file but need not ever reside in memory in that form.
c
        integer :: openfiles = -1 ! the number of open files
        integer :: h5err = 0      ! error status
                                  !  0 = OK
                                  ! -1 = Error
        integer(hid_t) :: H5T_NATIVE_DCMPLX  ! Data type for double
                                             ! complex entities.
        integer(hid_t) :: H5T_NATIVE_LOGICAL ! Data type for logical
                                             ! entities.
c
        integer(hid_t) :: proplistid ! property list identifier
c
        interface fwrite
          module procedure 
     &      fwrite_d0, fwrite_i0, fwrite_l0, fwrite_t0, fwrite_z0,
     &      fwrite_d1, fwrite_i1, fwrite_l1, fwrite_t1, fwrite_z1,
     &      fwrite_d2, fwrite_i2, fwrite_l2, fwrite_t2, fwrite_z2,
     &      fwrite_d3, fwrite_i3,            fwrite_t3, fwrite_z3,
     &      fwrite_d4, fwrite_i4,                       fwrite_z4,
     &      fwrite_d5,
     &      fwrite_d6,
     &      fwrite_d7
        end interface fwrite
c
        interface fread
          module procedure
     &      fread_d0,  fread_i0,  fread_l0,  fread_t0,  fread_z0,
     &      fread_d1,  fread_i1,  fread_l1,  fread_t1,  fread_z1,
     &      fread_d2,  fread_i2,  fread_l2,  fread_t2,  fread_z2,
     &      fread_d3,  fread_i3,             fread_t3,  fread_z3,
     &      fread_d4,  fread_i4,                        fread_z4,
     &      fread_d5,
     &      fread_d6,
     &      fread_d7
        end interface fread
c
      contains
c
c       ----------------------------------------------------------------
c
c       HDF5 provides a Fortran interface that does not support a
c       DOUBLE COMPLEX data type!? Hence this subroutine fixes that
c       omission.
c
        subroutine create_dcmplx_type()
        use hdf5
        implicit none
        integer(hsize_t) :: offset = 0
        integer          :: size8  = selected_real_kind(15)
        integer(hsize_t) :: size16 = 16
        call h5tcreate_f(H5T_COMPOUND_F,size16,H5T_NATIVE_DCMPLX,h5err)
        if (h5err.ne.0) then
          write(*,*)'ERROR: create_dcmplx_type: h5tcreate_f failed'
        endif
        call h5tinsert_f(H5T_NATIVE_DCMPLX,"real",offset,
     &                   h5kind_to_type(size8,H5_REAL_KIND),h5err)
        if (h5err.ne.0) then
          write(*,*)'ERROR: create_dcmplx_type: h5tinsert_f real failed'
        endif
        offset = offset + size8
        call h5tinsert_f(H5T_NATIVE_DCMPLX,"imag",offset,
     &                   h5kind_to_type(size8,h5_real_kind),h5err)
        if (h5err.ne.0) then
          write(*,*)'ERROR: create_dcmplx_type: h5tinsert_f imag failed'
        endif
        end subroutine create_dcmplx_type
c
c       ----------------------------------------------------------------
c
c       HDF5 provides a Fortran interface that does not support a
c       LOGICAL data type!? Hence this subroutine fixes that omission.
c
        subroutine create_logical_type()
        use hdf5
        implicit none
        integer          :: idummy
        integer(hsize_t) :: offset = 0
        integer(hsize_t) :: size4  = sizeof(idummy)
        call h5tcreate_f(H5T_OPAQUE_F,size4,H5T_NATIVE_LOGICAL,h5err)
        if (h5err.ne.0) then
          write(*,*)'ERROR: create_logical_type: h5tcreate_f failed'
        endif
        end subroutine create_logical_type
c
c       ----------------------------------------------------------------
c
c       Compute the number of elements in an N-dimensional array.
c
        integer function arraysize(nd,n)
        implicit none
        integer, intent(in) :: nd    ! the number of dimensions
        integer, intent(in) :: n(nd) ! the actual dimensions
        integer :: nsize
        integer :: ii ! counter
        nsize = n(1)
        do ii = 2, nd
          nsize = nsize * n(ii)
        enddo
        arraysize = nsize
        return
        end function arraysize
c
c       ----------------------------------------------------------------
c
c       Initialize the HDF5 library and set up the data types we need.
c
        subroutine finit()
        use hdf5
        implicit none
        if (openfiles.lt.0) then
          call h5open_f(h5err)
          if (h5err.ne.0) then
            write(*,*)'ERROR: finit: h5open_f failed'
          endif
          call create_dcmplx_type()
          call create_logical_type()
          openfiles = 0
        endif
        end subroutine finit
c
c       ----------------------------------------------------------------
c
c       Open a file for access. If no files are open yet then we
c       need to initialize HDF5 as well with a h5open_f call.
c
        subroutine fopen(filename,fileid)
        use hdf5
        implicit none
        character*(*), intent(in) :: filename 
        integer(hid_t), intent(out) :: fileid
        logical :: file_exists
        if (openfiles.eq.0) then
          ! create property list
          call h5pcreate_f(H5P_FILE_ACCESS_F,proplistid,h5err)
          if (h5err.ne.0) then
            write(*,*)'ERROR: fopen: h5pcreate_f failed'
          endif
          call h5pset_fapl_mpio_f(proplistid,MPI_COMM_WORLD,
     &                            MPI_INFO_NULL,h5err)
          if (h5err.ne.0) then
            write(*,*)'ERROR: fopen: h5pset_fapl_mpio_f failed'
          endif
        endif
        inquire(file=filename,exist=file_exists)
        if (file_exists) then
          call h5fopen_f(filename,H5F_ACC_RDWR_F,fileid,h5err,
     &                   access_prp = proplistid)
          if (h5err.ne.0) then
            write(*,*)'ERROR: fopen: h5fopen_f failed',filename
          endif
        else
          call h5fcreate_f(filename,H5F_ACC_EXCL_F,fileid,h5err,
     &                     access_prp = proplistid)
          if (h5err.ne.0) then
            write(*,*)'ERROR: fopen: h5fcreate_f failed',filename
          endif
        endif
        openfiles = openfiles + 1
        end subroutine fopen
c
c       ----------------------------------------------------------------
c
c       Close a file. If no more files are open then close the property
c       list as well.
c
        subroutine fclose(fileid)
        use hdf5
        implicit none
        integer(hid_t), intent(in) :: fileid
        call h5fclose_f(fileid,h5err)
        if (h5err.ne.0) then
          write(*,*)'ERROR: fclose: h5fclose_f failed'
        endif
        openfiles = openfiles - 1
        if (openfiles.le.0) then
          call h5pclose_f(proplistid,h5err)
          if (h5err.ne.0) then
            write(*,*)'ERROR: fclose: h5pclose_f failed'
          endif
        endif
        end subroutine fclose
c
c       ----------------------------------------------------------------
c
c       Close HDF5 down.
c
        subroutine ffinal()
        use hdf5
        implicit none
        if (openfiles.le.0) then
          call h5close_f(h5err)
          if (h5err.ne.0) then
            write(*,*)'ERROR: ffinal: h5close_f failed'
          endif
          openfiles = -1
        endif
        end subroutine ffinal
c
c       ----------------------------------------------------------------
c
c       Write an N-dimensional double precision array
c
        subroutine fwrite_dn(comm,darray,nd,ng,nl,nv,path,fileid)
        use hdf5
        implicit none
        integer, intent(in) :: comm   ! communicator for proc subset
        integer, intent(in) :: nd     ! number of dimensions
        integer, intent(in) :: ng(nd) ! the global dimension sizes
        integer, intent(in) :: nl(nd) ! the local dimension sizes
        integer, intent(in) :: nv(nd) ! the vector of offsets between
                                      ! global and local
        integer(hid_t), intent(in) :: fileid ! the file identifier
c
        character(len=*), intent(in) :: path ! where to store
c
        double precision, intent(in) :: darray(*) ! the data
c
c       Local
c
        integer(hsize_t) :: h5dimsf(1:nd) ! file space dimensions
        integer(hsize_t) :: h5dimsm(1:nd) ! memory space dimensions
        integer(hsize_t) :: h5dimsv(1:nd) ! offset vector
        integer(hsize_t) :: h5dimss(1:nd) ! stride (memory/file ???)
        integer(hsize_t) :: h5dimsc(1:nd) ! count: #blocks in each dim
        integer(hid_t)   :: h5fspaceid    ! file space ID
        integer(hid_t)   :: h5mspaceid    ! memory space ID
        integer(hid_t)   :: h5dsetid      ! data set ID
        integer(hid_t)   :: h5dmemtype    ! data memory type
        integer          :: h5select      ! selection parameter
c
        integer          :: ii            ! counter
c
        integer :: comm_size ! the size of communicator comm
        integer :: comm_rank ! the rank in communicator comm
        integer :: wrld_size ! the size of communicator MPI_COMM_WORLD
        integer :: wrld_rank ! the rank in communicator MPI_COMM_WORLD
        integer :: mpierr    ! the MPI error flag
c
        logical :: dataset_exists   ! does the data exist already?
        logical :: global_size_zero ! is the global dataset size 0?
        logical :: local_size_zero  ! is the local block size 0?
c
c       Code
c
        call mpi_comm_size(MPI_COMM_WORLD,wrld_size,mpierr)
        call mpi_comm_size(comm,comm_size,mpierr)
        call mpi_comm_rank(MPI_COMM_WORLD,wrld_rank,mpierr)
        call mpi_comm_rank(comm,comm_rank,mpierr)
c
        global_size_zero = .false.
        local_size_zero  = .false.
        do ii = 1, nd
          global_size_zero = global_size_zero.or.(ng(ii).le.0)
          local_size_zero  = local_size_zero .or.(nl(ii).le.0)
        enddo
        if (global_size_zero) then
          ! This dataset holds no data on any processor so just do
          ! nothing and return
          return
        endif
        h5dimsf = ng
        h5dimsm = nl
        h5dimsv = nv
        h5dimss = 1
        h5dimsc = 1
        h5dmemtype = H5T_NATIVE_DOUBLE
        h5select   = H5S_SELECT_SET_F
        call h5screate_simple_f(nd,h5dimsm,h5mspaceid,h5err)
        if (h5err.ne.0) then
          write(*,*)'ERROR: fwrite_dn: h5screate_simple_f failed:',path
        endif
c
        dataset_exists = fexist(fileid,path)
        if (dataset_exists) then
          call h5dopen_f(fileid,path,h5dsetid,h5err)
          if (h5err.ne.0) then
            write(*,*)'ERROR: fwrite_dn: h5dopen_f failed:',path
          endif
        else
          call h5screate_simple_f(nd,h5dimsf,h5fspaceid,h5err)
          if (h5err.ne.0) then
            write(*,*)'ERROR: fwrite_dn: h5screate_simple_f failed:',
     &                path
          endif
c         call h5dcreate_f(fileid,path,h5dmemtype,h5fspaceid,
c    &                     h5dsetid,h5err)
          call fcreate(fileid,path,h5dmemtype,h5fspaceid,h5dsetid)
          call h5sclose_f(h5fspaceid,h5err)
          if (h5err.ne.0) then
            write(*,*)'ERROR: fwrite_dn: h5sclose_f failed:',path
          endif
        endif
c
        call h5dget_space_f(h5dsetid,h5fspaceid,h5err)
        if (h5err.ne.0) then
          write(*,*)'ERROR: fwrite_dn: h5dget_space_f failed:',path
        endif
        call h5sselect_hyperslab_f(h5fspaceid,h5select,h5dimsv,h5dimsc,
     &                             h5err,h5dimss,h5dimsm)
        if (h5err.ne.0) then
          write(*,*)'ERROR: fwrite_dn: h5sselect_hyperslab_f failed:',
     &              path
        endif
c
        if (comm_rank.eq.0.and..not.local_size_zero) then
          call h5dwrite_f(h5dsetid,h5dmemtype,darray,h5dimsf,
     &                    h5err,file_space_id=h5fspaceid,
     &                    mem_space_id=h5mspaceid)
          if (h5err.ne.0) then
            write(*,*)'ERROR: fwrite_dn: h5dwrite_f failed:',path
          endif
        else
          ! Either we are not on rank 0, or the local data block size
          ! is zero. In either case let another processor handle the
          ! actual I/O and do nothing here. Note that the creation of
          ! the dataset needs to happen on all processors anyway to
          ! ensure that data layout on the file is the same on all
          ! processors.
        endif
        call h5dclose_f(h5dsetid,h5err)
        if (h5err.ne.0) then
          write(*,*)'ERROR: fwrite_dn: h5dclose_f failed:',path
        endif
        call h5sclose_f(h5fspaceid,h5err)
        if (h5err.ne.0) then
          write(*,*)'ERROR: fwrite_dn: h5sclose_f h5fspaceid failed:',
     &              path
        endif
        call h5sclose_f(h5mspaceid,h5err)
        if (h5err.ne.0) then
          write(*,*)'ERROR: fwrite_dn: h5sclose_f h5mspaceid failed:',
     &              path
        endif
c
        call mpi_barrier(comm,mpierr)
c
        end subroutine fwrite_dn
c
c       ----------------------------------------------------------------
c
c       Write an N-dimensional integer array
c
        subroutine fwrite_in(comm,iarray,nd,ng,nl,nv,path,fileid)
        use hdf5
        implicit none
        integer, intent(in) :: comm   ! communicator for proc subset
        integer, intent(in) :: nd     ! number of dimensions
        integer, intent(in) :: ng(nd) ! the global dimension sizes
        integer, intent(in) :: nl(nd) ! the local dimension sizes
        integer, intent(in) :: nv(nd) ! the vector of offsets between
                                      ! global and local
        integer(hid_t), intent(in) :: fileid ! the file identifier
c
        character(len=*), intent(in) :: path ! where to store
c
        integer, intent(in) :: iarray(*) ! the data
c
c       Local
c
        integer(hsize_t) :: h5dimsf(1:nd) ! file space dimensions
        integer(hsize_t) :: h5dimsm(1:nd) ! memory space dimensions
        integer(hsize_t) :: h5dimsv(1:nd) ! offset vector
        integer(hsize_t) :: h5dimss(1:nd) ! stride (file)
        integer(hsize_t) :: h5dimsc(1:nd) ! count: #blocks in each dim
        integer(hid_t)   :: h5fspaceid    ! file space ID
        integer(hid_t)   :: h5mspaceid    ! memory space ID
        integer(hid_t)   :: h5dsetid      ! data set ID
        integer(hid_t)   :: h5dmemtype    ! data memory type
        integer          :: h5select      ! selection parameter
c
        integer          :: ii            ! counter
c
        integer :: wrld_size ! the size of communicator MPI_COMM_WORLD
        integer :: wrld_rank ! the rank in communicator MPI_COMM_WORLD
        integer :: comm_size ! the size of communicator comm
        integer :: comm_rank ! the rank in communicator comm
        integer :: mpierr    ! the MPI error flag
c
        logical :: dataset_exists   ! does the data exist already?
        logical :: global_size_zero ! is the global dataset size 0?
        logical :: local_size_zero  ! is the local block size 0?
c
c       Code
c
        call mpi_comm_size(MPI_COMM_WORLD,wrld_size,mpierr)
        call mpi_comm_size(comm,comm_size,mpierr)
        call mpi_comm_rank(MPI_COMM_WORLD,wrld_rank,mpierr)
        call mpi_comm_rank(comm,comm_rank,mpierr)
c
        global_size_zero = .false.
        local_size_zero  = .false.
        do ii = 1, nd
          global_size_zero = global_size_zero.or.(ng(ii).le.0)
          local_size_zero  = local_size_zero .or.(nl(ii).le.0)
        enddo
        if (global_size_zero) then
          ! This dataset holds no data on any processor so just do
          ! nothing and return
          return
        endif
        h5dimsf = ng
        h5dimsm = nl
        h5dimsv = nv
        h5dimss = 1
        h5dimsc = 1
        h5dmemtype = H5T_NATIVE_INTEGER
        h5select   = H5S_SELECT_SET_F
        call h5screate_simple_f(nd,h5dimsm,h5mspaceid,h5err)
        if (h5err.ne.0) then
          write(*,*)'ERROR: fwrite_in: h5screate_simple_f failed:',path
        endif
c
c       call h5lexists_f(fileid,path,dataset_exists,h5err)
        dataset_exists = fexist(fileid,path)
        if (dataset_exists) then
          call h5dopen_f(fileid,path,h5dsetid,h5err)
          if (h5err.ne.0) then
            write(*,*)'ERROR: fwrite_in: h5dopen_f failed:',path
          endif
        else
          call h5screate_simple_f(nd,h5dimsf,h5fspaceid,h5err)
          if (h5err.ne.0) then
            write(*,*)'ERROR: fwrite_in: h5screate_simple_f failed:',
     &                path
          endif
c         call h5dcreate_f(fileid,path,h5dmemtype,h5fspaceid,
c    &                     h5dsetid,h5err)
          call fcreate(fileid,path,h5dmemtype,h5fspaceid,h5dsetid)
          call h5sclose_f(h5fspaceid,h5err)
          if (h5err.ne.0) then
            write(*,*)'ERROR: fwrite_in: h5sclose_f failed:',path
          endif
        endif
c
        call h5dget_space_f(h5dsetid,h5fspaceid,h5err)
        if (h5err.ne.0) then
          write(*,*)'ERROR: fwrite_in: h5dget_space_f failed:',path
        endif
        call h5sselect_hyperslab_f(h5fspaceid,h5select,h5dimsv,h5dimsc,
     &                             h5err,h5dimss,h5dimsm)
        if (h5err.ne.0) then
          write(*,*)'ERROR: fwrite_in: h5sselect_hyperslab_f failed:',
     &              path
        endif
c
        if (comm_rank.eq.0.and..not.local_size_zero) then
          call h5dwrite_f(h5dsetid,h5dmemtype,iarray,h5dimsf,
     &                    h5err,file_space_id=h5fspaceid,
     &                    mem_space_id=h5mspaceid)
          if (h5err.ne.0) then
            write(*,*)'ERROR: fwrite_in: h5dwrite_f failed:',path
          endif
        else
          ! Either we are not on rank 0, or the local data block size
          ! is zero. In either case let another processor handle the
          ! actual I/O and do nothing here. Note that the creation of
          ! the dataset needs to happen on all processors anyway to
          ! ensure that data layout on the file is the same on all
          ! processors.
        endif
        call h5dclose_f(h5dsetid,h5err)
        if (h5err.ne.0) then
          write(*,*)'ERROR: fwrite_in: h5dclose_f failed:',path
        endif
        call h5sclose_f(h5fspaceid,h5err)
        if (h5err.ne.0) then
          write(*,*)'ERROR: fwrite_in: h5sclose_f h5fspaceid failed:',
     &              path
        endif
        call h5sclose_f(h5mspaceid,h5err)
        if (h5err.ne.0) then
          write(*,*)'ERROR: fwrite_in: h5sclose_f h5mspaceid failed:',
     &              path
        endif
c
        call mpi_barrier(comm,mpierr)
c
        end subroutine fwrite_in
c
c       ----------------------------------------------------------------
c
c       Write an N-dimensional logical array
c
        subroutine fwrite_ln(comm,larray,nd,ng,nl,nv,path,fileid)
        use hdf5
        implicit none
        integer, intent(in) :: comm   ! communicator for proc subset
        integer, intent(in) :: nd     ! number of dimensions
        integer, intent(in) :: ng(nd) ! the global dimension sizes
        integer, intent(in) :: nl(nd) ! the local dimension sizes
        integer, intent(in) :: nv(nd) ! the vector of offsets between
                                      ! global and local
        integer(hid_t), intent(in) :: fileid ! the file identifier
c
        character(len=*), intent(in) :: path ! where to store
c
        logical, intent(in), target :: larray(*) ! the data
c
c       Local
c
        integer(hsize_t) :: h5dimsf(1:nd) ! file space dimensions
        integer(hsize_t) :: h5dimsm(1:nd) ! memory space dimensions
        integer(hsize_t) :: h5dimsv(1:nd) ! offset vector
        integer(hsize_t) :: h5dimss(1:nd) ! stride (file)
        integer(hsize_t) :: h5dimsc(1:nd) ! count: #blocks in each dim
        integer(hid_t)   :: h5fspaceid    ! file space ID
        integer(hid_t)   :: h5mspaceid    ! memory space ID
        integer(hid_t)   :: h5dsetid      ! data set ID
        integer(hid_t)   :: h5dmemtype    ! data memory type
        integer          :: h5select      ! selection parameter
        type(c_ptr)      :: f_ptr
c
        integer          :: ii            ! counter
c
        integer :: wrld_size ! the size of MPI_COMM_WORLD
        integer :: wrld_rank ! the rank in MPI_COMM_WORLD
        integer :: comm_size ! the size of communicator comm
        integer :: comm_rank ! the rank in communicator comm
        integer :: mpierr    ! the MPI error flag
c
        logical :: dataset_exists   ! does the data exist already?
        logical :: global_size_zero ! is the global dataset size 0?
        logical :: local_size_zero  ! is the local block size 0?
c
c       Code
c
        call mpi_comm_size(MPI_COMM_WORLD,wrld_size,mpierr)
        call mpi_comm_size(comm,comm_size,mpierr)
        call mpi_comm_rank(MPI_COMM_WORLD,wrld_rank,mpierr)
        call mpi_comm_rank(comm,comm_rank,mpierr)
c
        global_size_zero = .false.
        local_size_zero  = .false.
        do ii = 1, nd
          global_size_zero = global_size_zero.or.(ng(ii).le.0)
          local_size_zero  = local_size_zero .or.(nl(ii).le.0)
        enddo
        if (global_size_zero) then
          ! This dataset holds no data on any processor so just do
          ! nothing and return
          return
        endif
        h5dimsf = ng
        h5dimsm = nl
        h5dimsv = nv
        h5dimss = 1
        h5dimsc = 1
        h5dmemtype = H5T_NATIVE_LOGICAL
        h5select   = H5S_SELECT_SET_F
        call h5screate_simple_f(nd,h5dimsm,h5mspaceid,h5err)
        if (h5err.ne.0) then
          write(*,*)'ERROR: fwrite_ln: h5screate_simple_f failed:',path
        endif
c
c       call h5lexists_f(fileid,path,dataset_exists,h5err)
        dataset_exists = fexist(fileid,path)
        if (dataset_exists) then
          call h5dopen_f(fileid,path,h5dsetid,h5err)
          if (h5err.ne.0) then
            write(*,*)'ERROR: fwrite_ln: h5dopen_f failed:',path
          endif
        else
          call h5screate_simple_f(nd,h5dimsf,h5fspaceid,h5err)
          if (h5err.ne.0) then
            write(*,*)'ERROR: fwrite_ln: h5screate_simple_f failed:',
     &                path
          endif
c         call h5dcreate_f(fileid,path,h5dmemtype,h5fspaceid,
c    &                     h5dsetid,h5err)
          call fcreate(fileid,path,h5dmemtype,h5fspaceid,h5dsetid)
          call h5sclose_f(h5fspaceid,h5err)
          if (h5err.ne.0) then
            write(*,*)'ERROR: fwrite_ln: h5sclose_f failed:',path
          endif
        endif
c
        call h5dget_space_f(h5dsetid,h5fspaceid,h5err)
        if (h5err.ne.0) then
          write(*,*)'ERROR: fwrite_ln: h5dget_space_f failed:',path
        endif
        call h5sselect_hyperslab_f(h5fspaceid,h5select,h5dimsv,h5dimsc,
     &                             h5err,h5dimss,h5dimsm)
        if (h5err.ne.0) then
          write(*,*)'ERROR: fwrite_ln: h5sselect_hyperslab_f failed:',
     &              path
        endif
c
        if (comm_rank.eq.0.and..not.local_size_zero) then
          f_ptr = c_loc(larray)
c         call h5dwrite_f(h5dsetid,h5dmemtype,f_ptr,h5dimsf,
c    &                    h5err,file_space_id=h5fspaceid,
c    &                    mem_space_id=h5mspaceid)
          call h5dwrite_f(h5dsetid,h5dmemtype,f_ptr,h5err)
          if (h5err.ne.0) then
            write(*,*)'ERROR: fwrite_ln: h5dwrite_f failed:',path
          endif
        else
          ! Either we are not on rank 0, or the local data block size
          ! is zero. In either case let another processor handle the
          ! actual I/O and do nothing here. Note that the creation of
          ! the dataset needs to happen on all processors anyway to
          ! ensure that data layout on the file is the same on all
          ! processors.
        endif
        call h5dclose_f(h5dsetid,h5err)
        if (h5err.ne.0) then
          write(*,*)'ERROR: fwrite_ln: h5dclose_f failed:',path
        endif
        call h5sclose_f(h5fspaceid,h5err)
        if (h5err.ne.0) then
          write(*,*)'ERROR: fwrite_ln: h5sclose_f h5fspaceid failed:',
     &              path
        endif
        call h5sclose_f(h5mspaceid,h5err)
        if (h5err.ne.0) then
          write(*,*)'ERROR: fwrite_ln: h5sclose_f h5mspaceid failed:',
     &              path
        endif
c
        call mpi_barrier(comm,mpierr)
c
        end subroutine fwrite_ln
c
c       ----------------------------------------------------------------
c
c       Write an N-dimensional text array
c
        subroutine fwrite_tn(comm,tarray,nd,ng,nl,nv,path,fileid)
        use hdf5
        implicit none
        integer, intent(in) :: comm   ! communicator for proc set
        integer, intent(in) :: nd     ! number of dimensions
        integer, intent(in) :: ng(nd) ! the global dimension sizes
        integer, intent(in) :: nl(nd) ! the local dimension sizes
        integer, intent(in) :: nv(nd) ! the vector of offsets between
                                      ! global and local
        integer(hid_t), intent(in) :: fileid ! the file identifier
c
        character(len=*), intent(in) :: path ! where to store
c
        character(len=*), intent(in), target :: tarray(*) ! the data
c
c       Local
c
        integer(hsize_t) :: h5dimsf(1:nd)
        integer(hsize_t) :: h5dimsm(1:nd) ! memory space dimensions
        integer(hsize_t) :: h5dimsv(1:nd) ! offset vector
        integer(hsize_t) :: h5dimss(1:nd) ! stride (memory/file ???)
        integer(hsize_t) :: h5dimsc(1:nd) ! count: #blocks in each dim
        integer(hsize_t) :: sizetxt
        integer(hid_t)   :: h5mtype       ! the memory type
        integer(hid_t)   :: h5ftype       ! the file type
        integer(hid_t)   :: h5fspaceid    ! file space ID
        integer(hid_t)   :: h5mspaceid    ! memory space ID
        integer(hid_t)   :: h5dspaceid
        integer(hid_t)   :: h5dsetid
        integer          :: h5select      ! selection parameter
        integer(hid_t)   :: stringid
        type(c_ptr)      :: f_ptr
c
        integer          :: ii            ! counter
c
        integer :: comm_size ! the size of communicator comm
        integer :: comm_rank ! the rank in communicator comm
        integer :: wrld_size ! the size of communicator MPI_COMM_WORLD
        integer :: wrld_rank ! the rank in communicator MPI_COMM_WORLD
        integer :: mpierr
c
        logical :: dataset_exists   ! does the data exist already?
        logical :: global_size_zero ! is the global dataset size 0?
        logical :: local_size_zero  ! is the local block size 0?
c
c       Code
c
        call mpi_comm_size(MPI_COMM_WORLD,wrld_size,mpierr)
        call mpi_comm_size(comm,comm_size,mpierr)
        call mpi_comm_rank(MPI_COMM_WORLD,wrld_rank,mpierr)
        call mpi_comm_rank(comm,comm_rank,mpierr)
c
        global_size_zero = .false.
        local_size_zero  = .false.
        do ii = 1, nd
          global_size_zero = global_size_zero.or.(ng(ii).le.0)
          local_size_zero  = local_size_zero .or.(nl(ii).le.0)
        enddo
        if (global_size_zero) then
          ! This dataset holds no data on any processor so just do
          ! nothing and return
          return
        endif
        h5dimsf = ng
        h5dimsm = nl
        h5dimsv = nv
        h5dimss = 1
        h5dimsc = 1
        sizetxt = len(tarray(1))
        stringid = H5T_FORTRAN_S1
        call h5tcopy_f(stringid,h5mtype,h5err)
        if (h5err.ne.0) then
          write(*,*)'ERROR: fwrite_tn: h5tcopy_f h5mtype failed:',path
        endif
        call h5tcopy_f(stringid,h5ftype,h5err)
        if (h5err.ne.0) then
          write(*,*)'ERROR: fwrite_tn: h5tcopy_f h5ftype failed:',path
        endif
        call h5tset_size_f(h5mtype,sizetxt,h5err)
        if (h5err.ne.0) then
          write(*,*)'ERROR: fwrite_tn: h5tset_size_f h5mtype failed:',
     &              path
        endif
        call h5tset_size_f(h5ftype,sizetxt,h5err)
        if (h5err.ne.0) then
          write(*,*)'ERROR: fwrite_tn: h5tset_size_f h5ftype failed:',
     &              path
        endif
        call h5screate_simple_f(nd,h5dimsm,h5dspaceid,h5err)
        if (h5err.ne.0) then
          write(*,*)'ERROR: fwrite_tn: h5screate_simple_f failed:',path
        endif
        dataset_exists = fexist(fileid,path)
        if (dataset_exists) then
          call h5dopen_f(fileid,path,h5dsetid,h5err)
          if (h5err.ne.0) then
            write(*,*)'ERROR: fwrite_tn: h5dopen_f failed:',path
          endif
        else
          call h5screate_simple_f(nd,h5dimsf,h5fspaceid,h5err)
          if (h5err.ne.0) then
            write(*,*)'ERROR: fwrite_tn: h5screate_simple_f failed:',
     &                path
          endif
c         call h5dcreate_f(fileid,path,stringid,h5dspaceid,
c    &                     h5dsetid,h5err)
          call fcreate(fileid,path,h5ftype,h5dspaceid,h5dsetid)
          call h5sclose_f(h5dspaceid,h5err)
          if (h5err.ne.0) then
            write(*,*)'ERROR: fwrite_tn: h5sclose_f failed:',path
          endif
        endif
        if (comm_rank.eq.0.and..not.local_size_zero) then
          f_ptr = c_loc(tarray)
          call h5dwrite_f(h5dsetid,h5mtype,f_ptr,h5err)
          if (h5err.ne.0) then
            write(*,*)'ERROR: fwrite_tn: h5dwrite_f failed:',path
          endif
        else
          ! Either we are not on rank 0, or the local data block size
          ! is zero. In either case let another processor handle the
          ! actual I/O and do nothing here. Note that the creation of
          ! the dataset needs to happen on all processors anyway to
          ! ensure that data layout on the file is the same on all
          ! processors.
        endif
        call h5dclose_f(h5dsetid,h5err)
        if (h5err.ne.0) then
          write(*,*)'ERROR: fwrite_tn: h5dclose_f failed:',path
        endif
        call h5tclose_f(h5mtype,h5err) 
        if (h5err.ne.0) then
          write(*,*)'ERROR: fwrite_tn: h5tclose_f h5mtype failed:',path
        endif
        call h5tclose_f(h5ftype,h5err) 
        if (h5err.ne.0) then
          write(*,*)'ERROR: fwrite_tn: h5tclose_f h5ftype failed:',path
        endif
c
        end subroutine fwrite_tn
c
c       ----------------------------------------------------------------
c
c       Write an N-dimensional double complex array
c
        subroutine fwrite_zn(comm,zarray,nd,ng,nl,nv,path,fileid)
        use hdf5
        implicit none
        integer, intent(in) :: comm   ! communicator for proc subset
        integer, intent(in) :: nd     ! number of dimensions
        integer, intent(in) :: ng(nd) ! the global dimension sizes
        integer, intent(in) :: nl(nd) ! the local dimension sizes
        integer, intent(in) :: nv(nd) ! the vector of offsets between
                                      ! global and local
        integer(hid_t), intent(in) :: fileid ! the file identifier
c
        character(len=*), intent(in) :: path ! where to store
c
        double complex, intent(in), target :: zarray(*) ! the data
c
c       Local
c
        integer(hsize_t) :: h5dimsf(1:nd) ! file space dimensions
        integer(hsize_t) :: h5dimsm(1:nd) ! memory space dimensions
        integer(hsize_t) :: h5dimsv(1:nd) ! offset vector
        integer(hsize_t) :: h5dimss(1:nd) ! stride (memory/file ???)
        integer(hsize_t) :: h5dimsc(1:nd) ! count: #blocks in each dim
        integer(hid_t)   :: h5fspaceid    ! file space ID
        integer(hid_t)   :: h5mspaceid    ! memory space ID
        integer(hid_t)   :: h5dsetid      ! data set ID
        integer(hid_t)   :: h5dmemtype    ! data memory type
        integer          :: h5select      ! selection parameter
c
        integer          :: ii            ! counter
c
        integer :: comm_size ! the size of communicator comm
        integer :: comm_rank ! the rank in communicator comm
        integer :: wrld_size ! the size of communicator MPI_COMM_WORLD
        integer :: wrld_rank ! the rank in communicator MPI_COMM_WORLD
        integer :: mpierr    ! the MPI error flag
c
        logical :: dataset_exists   ! does the data exist already?
        logical :: global_size_zero ! is the global dataset size 0?
        logical :: local_size_zero  ! is the local block size 0?
c
        type(c_ptr)      :: f_ptr
c
c       Code
c
c       Look at hdf5-1.8.19/fortran/examples/compound_complex_fortran2003.f90
c       for an example as to how to do this.
c
        call mpi_comm_size(MPI_COMM_WORLD,wrld_size,mpierr)
        call mpi_comm_size(comm,comm_size,mpierr)
        call mpi_comm_rank(MPI_COMM_WORLD,wrld_rank,mpierr)
        call mpi_comm_rank(comm,comm_rank,mpierr)
c
        global_size_zero = .false.
        local_size_zero  = .false.
        do ii = 1, nd
          global_size_zero = global_size_zero.or.(ng(ii).le.0)
          local_size_zero  = local_size_zero .or.(nl(ii).le.0)
        enddo
        if (global_size_zero) then
          ! This dataset holds no data on any processor so just do
          ! nothing and return
          return
        endif
        h5dimsf = ng
        h5dimsm = nl
        h5dimsv = nv
        h5dimss = 1
        h5dimsc = 1
        h5dmemtype = H5T_NATIVE_DCMPLX
        h5select   = H5S_SELECT_SET_F
        call h5screate_simple_f(nd,h5dimsm,h5mspaceid,h5err)
        if (h5err.ne.0) then
          write(*,*)'ERROR: fwrite_zn: h5screate_simple_f failed:',path
        endif
c
c       call h5lexists_f(fileid,path,dataset_exists,h5err)
        dataset_exists = fexist(fileid,path)
        if (dataset_exists) then
          call h5dopen_f(fileid,path,h5dsetid,h5err)
          if (h5err.ne.0) then
            write(*,*)'ERROR: fwrite_zn: h5dopen_f failed:',path
          endif
        else
          call h5screate_simple_f(nd,h5dimsf,h5fspaceid,h5err)
          if (h5err.ne.0) then
            write(*,*)'ERROR: fwrite_zn: h5screate_simple_f failed:',
     &                path
          endif
c         call h5dcreate_f(fileid,path,h5dmemtype,h5fspaceid,
c    &                     h5dsetid,h5err)
          call fcreate(fileid,path,h5dmemtype,h5fspaceid,h5dsetid)
          call h5sclose_f(h5fspaceid,h5err)
          if (h5err.ne.0) then
            write(*,*)'ERROR: fwrite_zn: h5sclose_f failed:',path
          endif
        endif
c
        call h5dget_space_f(h5dsetid,h5fspaceid,h5err)
        if (h5err.ne.0) then
          write(*,*)'ERROR: fwrite_zn: h5dget_space_f failed:',path
        endif
        call h5sselect_hyperslab_f(h5fspaceid,h5select,h5dimsv,h5dimsc,
     &                             h5err,h5dimss,h5dimsm)
        if (h5err.ne.0) then
          write(*,*)'ERROR: fwrite_zn: h5sselect_hyperslab_f failed:',
     &              path
        endif
c
        if (comm_rank.eq.0.and..not.local_size_zero) then
          f_ptr = c_loc(zarray)
          call h5dwrite_f(h5dsetid,h5dmemtype,f_ptr,
     &                    h5err,mem_space_id=h5mspaceid,
     &                    file_space_id=h5fspaceid)
          if (h5err.ne.0) then
            write(*,*)'ERROR: fwrite_zn: h5dwrite_f failed:',path
          endif
        else
          ! Either we are not on rank 0, or the local data block size
          ! is zero. In either case let another processor handle the
          ! actual I/O and do nothing here. Note that the creation of
          ! the dataset needs to happen on all processors anyway to
          ! ensure that data layout on the file is the same on all
          ! processors.
        endif
        call h5dclose_f(h5dsetid,h5err)
        if (h5err.ne.0) then
          write(*,*)'ERROR: fwrite_zn: h5dclose_f failed:',path
        endif
        call h5sclose_f(h5fspaceid,h5err)
        if (h5err.ne.0) then
          write(*,*)'ERROR: fwrite_zn: h5sclose_f h5fspaceid failed:',
     &              path
        endif
        call h5sclose_f(h5mspaceid,h5err)
        if (h5err.ne.0) then
          write(*,*)'ERROR: fwrite_zn: h5sclose_f h5mspaceid failed:',
     &              path
        endif
c
        call mpi_barrier(comm,mpierr)
c
        end subroutine fwrite_zn
c
c       ----------------------------------------------------------------
c
c       Read an N-dimensional double precision array
c
        subroutine fread_dn(comm,f_ptr,nd,ng,nl,nv,path,fileid)
        use hdf5
        implicit none
        integer, intent(in) :: comm   ! communicator for proc subset
        integer, intent(in) :: nd     ! number of dimensions
        integer, intent(in) :: ng(nd) ! the global dimension sizes
        integer, intent(in) :: nl(nd) ! the local dimension sizes
        integer, intent(in) :: nv(nd) ! the vector of offsets between
                                      ! global and local
        integer(hid_t), intent(in) :: fileid ! the file identifier
c       
        character(len=*), intent(in) :: path ! where to store
c
        type(c_ptr) :: f_ptr ! pointer to the data
c
c       Local
c
        integer(hsize_t) :: h5dimsf(1:nd) ! file space dimensions
        integer(hsize_t) :: h5dimsm(1:nd) ! memory space dimensions
        integer(hsize_t) :: h5dimsv(1:nd) ! offset vector
        integer(hsize_t) :: h5dimss(1:nd) ! stride (file)
        integer(hsize_t) :: h5dimsc(1:nd) ! count: #blocks in each dim
        integer(hid_t)   :: h5fspaceid    ! file space ID
        integer(hid_t)   :: h5mspaceid    ! memory space ID
        integer(hid_t)   :: h5dsetid      ! data set ID
        integer(hid_t)   :: h5dmemtype    ! data memory type
        integer          :: h5select      ! selection parameter
c
        integer          :: ii            ! counter
c
        integer :: comm_size ! the size of communicator comm
        integer :: comm_rank ! the rank in communicator comm
        integer :: wrld_size ! the size of communicator MPI_COMM_WORLD
        integer :: wrld_rank ! the rank in communicator MPI_COMM_WORLD
        integer :: mpierr    ! the MPI error flag
c
        logical :: dataset_exists   ! does the data exist already?
        logical :: global_size_zero ! is the global dataset size 0?
        logical :: local_size_zero  ! is the local block size 0?
c
c       Code
c
        call mpi_comm_size(MPI_COMM_WORLD,wrld_size,mpierr)
        call mpi_comm_size(comm,comm_size,mpierr)
        call mpi_comm_rank(MPI_COMM_WORLD,wrld_rank,mpierr)
        call mpi_comm_rank(comm,comm_rank,mpierr)
c
        global_size_zero = .false.
        local_size_zero  = .false.
        do ii = 1, nd
          global_size_zero = global_size_zero.or.(ng(ii).le.0)
          local_size_zero  = local_size_zero .or.(nl(ii).le.0)
        enddo
        if (global_size_zero) then
          ! This dataset holds no data on any processor so just do
          ! nothing and return
          return
        endif
        h5dimsf = ng
        h5dimsm = nl
        h5dimsv = nv
        h5dimss = 1
        h5dimsc = 1
        h5dmemtype = H5T_NATIVE_DOUBLE
        h5select   = H5S_SELECT_SET_F
        call h5screate_simple_f(nd,h5dimsm,h5mspaceid,h5err)
        if (h5err.ne.0) then
          write(*,*)'ERROR: fread_dn: h5screate_simple_f failed',path
        endif
c
c       call h5lexists_f(fileid,path,dataset_exists,h5err)
        dataset_exists = fexist(fileid,path)
        if (.not.dataset_exists) then
          write(*,*)wrld_rank," dataset: ",path," does not exist!"
          flush(6)
          call h5sclose_f(h5mspaceid,h5err)
          if (h5err.ne.0) then
            write(*,*)'ERROR: fread_dn: h5sclose_f failed',path
          endif
          return
        endif
        call h5dopen_f(fileid,path,h5dsetid,h5err)
        if (h5err.ne.0) then
          write(*,*)'ERROR: fread_dn: h5dopen_f failed',path
        endif
c
        call h5dget_space_f(h5dsetid,h5fspaceid,h5err)
        if (h5err.ne.0) then
          write(*,*)'ERROR: fread_dn: h5dget_space_f failed',path
        endif
        call h5sselect_hyperslab_f(h5fspaceid,h5select,h5dimsv,h5dimsc,
     &                             h5err,h5dimss,h5dimsm)
        if (h5err.ne.0) then
          write(*,*)'ERROR: fread_dn: h5sselect_hyperslab_f failed',path
        endif
c
        if (.not.local_size_zero) then
          call h5dread_f(h5dsetid,h5dmemtype,f_ptr,
     &                   h5err,mem_space_id=h5mspaceid,
     &                   file_space_id=h5fspaceid)
          if (h5err.ne.0) then
            write(*,*)'ERROR: fread_dn: h5dread_f failed',path
          endif
        else
          ! The local block size is zero, so no I/O can be done here.
          ! The other operations on the file check the consistency
          ! of the file across a number of processors. The data layout
          ! should be the same on all processors.
        endif
c
        if (dataset_exists) then
          call h5dclose_f(h5dsetid,h5err)
          if (h5err.ne.0) then
            write(*,*)'ERROR: fread_dn: h5dclose_f failed',path
          endif
        endif
        call h5sclose_f(h5fspaceid,h5err)
        if (h5err.ne.0) then
          write(*,*)'ERROR: fread_dn: h5sclose_f h5fspaceid failed',path
        endif
        call h5sclose_f(h5mspaceid,h5err)
        if (h5err.ne.0) then
          write(*,*)'ERROR: fread_dn: h5sclose_f h5mspaceid failed',path
        endif
c
        call mpi_barrier(comm,mpierr)
c
        end subroutine fread_dn
c
c       ----------------------------------------------------------------
c
c       Read an N-dimensional integer array
c
        subroutine fread_in(comm,f_ptr,nd,ng,nl,nv,path,fileid)
        use hdf5
        implicit none
        integer, intent(in) :: comm   ! communicator for proc subset
        integer, intent(in) :: nd     ! number of dimensions
        integer, intent(in) :: ng(nd) ! the global dimension sizes
        integer, intent(in) :: nl(nd) ! the local dimension sizes
        integer, intent(in) :: nv(nd) ! the vector of offsets between
                                      ! global and local
        integer(hid_t), intent(in) :: fileid ! the file identifier
c       
        character(len=*), intent(in) :: path ! where to store
c
        type(c_ptr) :: f_ptr ! the location of the data
c
c       Local
c
        integer(hsize_t) :: h5dimsf(1:nd) ! file space dimensions
        integer(hsize_t) :: h5dimsm(1:nd) ! memory space dimensions
        integer(hsize_t) :: h5dimsv(1:nd) ! offset vector
        integer(hsize_t) :: h5dimss(1:nd) ! stride (file)
        integer(hsize_t) :: h5dimsc(1:nd) ! count: #blocks in each dim
        integer(hid_t)   :: h5fspaceid    ! file space ID
        integer(hid_t)   :: h5mspaceid    ! memory space ID
        integer(hid_t)   :: h5dsetid      ! data set ID
        integer(hid_t)   :: h5dmemtype    ! data memory type
        integer          :: h5select      ! selection parameter
c
        integer          :: ii            ! counter
c
        integer :: wrld_rank ! the rank in MPI_COMM_WORLD
        integer :: mpierr    ! the MPI error flag
c
        logical :: dataset_exists   ! does the data exist already?
        logical :: global_size_zero ! is the global dataset size 0?
        logical :: local_size_zero  ! is the local block size 0?
c
c       Code
c
        call mpi_comm_rank(MPI_COMM_WORLD,wrld_rank,mpierr)
c
        global_size_zero = .false.
        local_size_zero  = .false.
        do ii = 1, nd
          global_size_zero = global_size_zero.or.(ng(ii).le.0)
          local_size_zero  = local_size_zero .or.(nl(ii).le.0)
        enddo
        if (global_size_zero) then
          ! This dataset holds no data on any processor so just do
          ! nothing and return
          return
        endif
        h5dimsf = ng
        h5dimsm = nl
        h5dimsv = nv
        h5dimss = 1
        h5dimsc = 1
        h5dmemtype = H5T_NATIVE_INTEGER
        h5select   = H5S_SELECT_SET_F
        call h5screate_simple_f(nd,h5dimsm,h5mspaceid,h5err)
        if (h5err.ne.0) then
          write(*,*)'ERROR: fread_in: h5screate_simple_f failed',path
        endif
c
c       call h5lexists_f(fileid,path,dataset_exists,h5err)
        dataset_exists = fexist(fileid,path)
        if (.not.dataset_exists) then
          write(*,*)wrld_rank," dataset: ",path," does not exist!"
          call h5sclose_f(h5mspaceid,h5err)
          if (h5err.ne.0) then
            write(*,*)'ERROR: fread_in: h5sclose_f failed',path
          endif
          return
        endif
        call h5dopen_f(fileid,path,h5dsetid,h5err)
        if (h5err.ne.0) then
          write(*,*)'ERROR: fread_in: h5dopen_f failed',path
        endif
c
        call h5dget_space_f(h5dsetid,h5fspaceid,h5err)
        if (h5err.ne.0) then
          write(*,*)'ERROR: fread_in: h5dget_space_f failed',path
        endif
        call h5sselect_hyperslab_f(h5fspaceid,h5select,h5dimsv,h5dimsc,
     &                             h5err,h5dimss,h5dimsm)
        if (h5err.ne.0) then
          write(*,*)'ERROR: fread_in: h5sselect_hyperslab_f failed',path
        endif
c
        if (.not.local_size_zero) then
          call h5dread_f(h5dsetid,h5dmemtype,f_ptr,
     &                   h5err,mem_space_id=h5mspaceid,
     &                   file_space_id=h5fspaceid)
          if (h5err.ne.0) then
            write(*,*)'ERROR: fread_in: h5dread_f failed',path
          endif
        else
          ! The local block size is zero, so no I/O can be done here.
          ! The other operations on the file check the consistency
          ! of the file across a number of processors. The data layout
          ! should be the same on all processors.
        endif
c
        if (dataset_exists) then
          call h5dclose_f(h5dsetid,h5err)
          if (h5err.ne.0) then
            write(*,*)'ERROR: fread_in: h5dclose_f failed',path
          endif
        endif
        call h5sclose_f(h5fspaceid,h5err)
        if (h5err.ne.0) then
          write(*,*)'ERROR: fread_in: h5sclose_f h5fspaceid failed',path
        endif
        call h5sclose_f(h5mspaceid,h5err)
        if (h5err.ne.0) then
          write(*,*)'ERROR: fread_in: h5sclose_f h5mspaceid failed',path
        endif
c
        call mpi_barrier(comm,mpierr)
c
        end subroutine fread_in
c
c       ----------------------------------------------------------------
c
c       Read an N-dimensional logical array
c
        subroutine fread_ln(comm,f_ptr,nd,ng,nl,nv,path,fileid)
        use hdf5
        implicit none
        integer, intent(in) :: comm   ! communicator for proc subset
        integer, intent(in) :: nd     ! number of dimensions
        integer, intent(in) :: ng(nd) ! the global dimension sizes
        integer, intent(in) :: nl(nd) ! the local dimension sizes
        integer, intent(in) :: nv(nd) ! the vector of offsets between
                                      ! global and local
        integer(hid_t), intent(in) :: fileid ! the file identifier
c       
        character(len=*), intent(in) :: path ! where to store
c
        type(c_ptr) :: f_ptr ! the location of the data
c
c       Local
c
        integer(hsize_t) :: h5dimsf(1:nd) ! file space dimensions
        integer(hsize_t) :: h5dimsm(1:nd) ! memory space dimensions
        integer(hsize_t) :: h5dimsv(1:nd) ! offset vector
        integer(hsize_t) :: h5dimss(1:nd) ! stride (file)
        integer(hsize_t) :: h5dimsc(1:nd) ! count: #blocks in each dim
        integer(hid_t)   :: h5fspaceid    ! file space ID
        integer(hid_t)   :: h5mspaceid    ! memory space ID
        integer(hid_t)   :: h5dsetid      ! data set ID
        integer(hid_t)   :: h5dmemtype    ! data memory type
        integer          :: h5select      ! selection parameter
c
        integer          :: ii            ! counter
c
        integer :: wrld_rank ! the rank in MPI_COMM_WORLD
        integer :: mpierr    ! the MPI error flag
c
        logical :: dataset_exists   ! does the data exist already?
        logical :: global_size_zero ! is the global dataset size 0?
        logical :: local_size_zero  ! is the local block size 0?
c
c       Code
c
        call mpi_comm_rank(MPI_COMM_WORLD,wrld_rank,mpierr)
c
        global_size_zero = .false.
        local_size_zero  = .false.
        do ii = 1, nd
          global_size_zero = global_size_zero.or.(ng(ii).le.0)
          local_size_zero  = local_size_zero .or.(nl(ii).le.0)
        enddo
        if (global_size_zero) then
          ! This dataset holds no data on any processor so just do
          ! nothing and return
          return
        endif
        h5dimsf = ng
        h5dimsm = nl
        h5dimsv = nv
        h5dimss = 1
        h5dimsc = 1
        h5dmemtype = H5T_NATIVE_LOGICAL
        h5select   = H5S_SELECT_SET_F
        call h5screate_simple_f(nd,h5dimsm,h5mspaceid,h5err)
        if (h5err.ne.0) then
          write(*,*)'ERROR: fread_ln: h5screate_simple_f failed',path
        endif
c
c       call h5lexists_f(fileid,path,dataset_exists,h5err)
        dataset_exists = fexist(fileid,path)
        if (.not.dataset_exists) then
          write(*,*)wrld_rank," dataset: ",path," does not exist!"
          call h5sclose_f(h5mspaceid,h5err)
          if (h5err.ne.0) then
            write(*,*)'ERROR: fread_ln: h5sclose_f failed',path
          endif
          return
        endif
        call h5dopen_f(fileid,path,h5dsetid,h5err)
        if (h5err.ne.0) then
          write(*,*)'ERROR: fread_ln: h5dopen_f failed',path
        endif
c
        call h5dget_space_f(h5dsetid,h5fspaceid,h5err)
        if (h5err.ne.0) then
          write(*,*)'ERROR: fread_ln: h5dget_space_f failed',path
        endif
        call h5sselect_hyperslab_f(h5fspaceid,h5select,h5dimsv,h5dimsc,
     &                             h5err,h5dimss,h5dimsm)
        if (h5err.ne.0) then
          write(*,*)'ERROR: fread_ln: h5sselect_hyperslab_f failed',path
        endif
c
        if (.not.local_size_zero) then
          call h5dread_f(h5dsetid,h5dmemtype,f_ptr,
     &                   h5err,mem_space_id=h5mspaceid,
     &                   file_space_id=h5fspaceid)
          if (h5err.ne.0) then
            write(*,*)'ERROR: fread_ln: h5dread_f failed',path
          endif
        else
          ! The local block size is zero, so no I/O can be done here.
          ! The other operations on the file check the consistency
          ! of the file across a number of processors. The data layout
          ! should be the same on all processors.
        endif
c
        if (dataset_exists) then
          call h5dclose_f(h5dsetid,h5err)
          if (h5err.ne.0) then
            write(*,*)'ERROR: fread_ln: h5dclose_f failed',path
          endif
        endif
        call h5sclose_f(h5fspaceid,h5err)
        if (h5err.ne.0) then
          write(*,*)'ERROR: fread_ln: h5sclose_f h5fspaceid failed',path
        endif
        call h5sclose_f(h5mspaceid,h5err)
        if (h5err.ne.0) then
          write(*,*)'ERROR: fread_ln: h5mclose_f h5fspaceid failed',path
        endif
c
        call mpi_barrier(comm,mpierr)
c
        end subroutine fread_ln
c
c       ----------------------------------------------------------------
c
c       Read an N-dimensional text array
c
        subroutine fread_tn(comm,f_ptr,nt,nd,ng,nl,nv,path,fileid)
        use hdf5
        implicit none
        integer, intent(in) :: comm ! communicator for proc set
        integer, intent(in) :: nt     ! number of chars in string
        integer, intent(in) :: nd     ! number of dimensions
        integer, intent(in) :: ng(nd) ! the global dimension sizes
        integer, intent(in) :: nl(nd) ! the local dimension sizes
        integer, intent(in) :: nv(nd) ! the vector of offsets between
                                      ! global and local
        integer(hid_t), intent(in) :: fileid ! the file identifier
c
        character(len=*), intent(in) :: path ! where to store
c
        type(c_ptr) :: f_ptr ! pointer to the data
c
c       Local
c
        integer(hsize_t) :: h5dimsf(1:nd) ! file space dimensions
        integer(hsize_t) :: h5dimsm(1:nd) ! memory space dimensions
        integer(hsize_t) :: h5dimsv(1:nd) ! offset vector
        integer(hsize_t) :: h5dimss(1:nd) ! stride (file)
        integer(hsize_t) :: h5dimsc(1:nd) ! count: #blocks in each dim
        integer(hsize_t) :: sizetxt       ! size of a single string
        integer(hid_t)   :: h5mtype       ! memory type
        integer(hid_t)   :: h5ftype       ! file type
        integer(hid_t)   :: h5fspaceid    ! file space ID
        integer(hid_t)   :: h5mspaceid    ! memory space ID
        integer(hid_t)   :: h5dspaceid
        integer(hid_t)   :: h5dsetid
        integer(hid_t)   :: stringid
        integer          :: h5select      ! selection parameter
c
        integer          :: ii            ! counter
c
        integer :: comm_size ! the size of communicator comm
        integer :: comm_rank ! the rank in communicator comm
        integer :: wrld_size ! the size of communicator MPI_COMM_WORLD
        integer :: wrld_rank ! the rank in communicator MPI_COMM_WORLD
        integer :: mpierr
c
        logical :: dataset_exists   ! does the data exist already?
        logical :: global_size_zero ! is the global dataset size 0?
        logical :: local_size_zero  ! is the local block size 0?
c
c       Code
c
        call mpi_comm_size(MPI_COMM_WORLD,wrld_size,mpierr)
        call mpi_comm_size(comm,comm_size,mpierr)
        call mpi_comm_rank(MPI_COMM_WORLD,wrld_rank,mpierr)
        call mpi_comm_rank(comm,comm_rank,mpierr)
c
        global_size_zero = .false.
        local_size_zero  = .false.
        do ii = 1, nd
          global_size_zero = global_size_zero.or.(ng(ii).le.0)
          local_size_zero  = local_size_zero .or.(nl(ii).le.0)
        enddo
        if (global_size_zero) then
          ! This dataset holds no data on any processor so just do
          ! nothing and return
          return
        endif
        sizetxt = nt
        h5dimsf = ng
        h5dimsm = nl
        h5dimsv = nv
        h5dimss = 1
        h5dimsc = 1
        h5select = H5S_SELECT_SET_F
        stringid = H5T_FORTRAN_S1
        call h5tcopy_f(stringid,h5mtype,h5err)
        if (h5err.ne.0) then
          write(*,*)'ERROR: fread_tn: h5tcopy_f failed',path
        endif
        call h5tset_size_f(h5mtype,sizetxt,h5err)
        if (h5err.ne.0) then
          write(*,*)'ERROR: fread_tn: h5tset_size_f failed',path
        endif
c
        call h5screate_simple_f(nd,h5dimsm,h5mspaceid,h5err)
        if (h5err.ne.0) then
          write(*,*)'ERROR: fread_tn: h5screate_simple_f failed',path
        endif
c
c       call h5lexists_f(fileid,path,dataset_exists,h5err)
        dataset_exists = fexist(fileid,path)
        if (.not.dataset_exists) then
          write(*,*)wrld_rank," dataset: ",path," does not exist!"
c         call h5sclose_f(h5mspaceid,h5err)
          return
        endif
        call h5dopen_f(fileid,path,h5dsetid,h5err)
        if (h5err.ne.0) then
          write(*,*)'ERROR: fread_tn: h5dopen_f failed',path
        endif
c
        call h5dget_space_f(h5dsetid,h5fspaceid,h5err)
        if (h5err.ne.0) then
          write(*,*)'ERROR: fread_tn: h5dget_space_f failed',path
        endif
        call h5sselect_hyperslab_f(h5fspaceid,h5select,h5dimsv,h5dimsc,
     &                             h5err,h5dimss,h5dimsm)
        if (h5err.ne.0) then
          write(*,*)'ERROR: fread_tn: h5sselect_hyperslab_f failed',path
        endif
c       call h5dget_type_f(h5dsetid,h5ftype,h5err)
c       call h5dget_space_f(h5dsetid,h5dspaceid,h5err)
c
        if (.not.local_size_zero) then
          call h5dread_f(h5dsetid,h5mtype,f_ptr,
     &                   h5err,mem_space_id=h5mspaceid,
     &                   file_space_id=h5fspaceid)
          if (h5err.ne.0) then
            write(*,*)'ERROR: fread_tn: h5dread_f failed',path
          endif
        else
          ! The local block size is zero, so no I/O can be done here.
          ! The other operations on the file check the consistency
          ! of the file across a number of processors. The data layout
          ! should be the same on all processors.
        endif
c
        if (dataset_exists) then
          call h5dclose_f(h5dsetid,h5err)
          if (h5err.ne.0) then
            write(*,*)'ERROR: fread_tn: h5dclose_f failed',path
          endif
        endif
        call h5sclose_f(h5fspaceid,h5err)
        if (h5err.ne.0) then
          write(*,*)'ERROR: fread_tn: h5sclose_f h5fspaceid failed',path
        endif
        call h5sclose_f(h5mspaceid,h5err)
        if (h5err.ne.0) then
          write(*,*)'ERROR: fread_tn: h5sclose_f h5mspaceid failed',path
        endif
c       call h5sclose_f(h5dspaceid,h5err)
c       call h5tclose_f(h5ftype,h5err)
        call h5tclose_f(h5mtype,h5err)
        if (h5err.ne.0) then
          write(*,*)'ERROR: fread_tn: h5tclose_f failed',path
        endif
c
        end subroutine fread_tn
c
c       ----------------------------------------------------------------
c
c       Read an N-dimensional double complex array
c
        subroutine fread_zn(comm,f_ptr,nd,ng,nl,nv,path,fileid)
        use hdf5
        implicit none
        integer, intent(in) :: comm   ! communicator for proc subset
        integer, intent(in) :: nd     ! number of dimensions
        integer, intent(in) :: ng(nd) ! the global dimension sizes
        integer, intent(in) :: nl(nd) ! the local dimension sizes
        integer, intent(in) :: nv(nd) ! the vector of offsets between
                                      ! global and local
        integer(hid_t), intent(in) :: fileid ! the file identifier
c       
        character(len=*), intent(in) :: path ! where to store
c
        type(c_ptr) :: f_ptr ! the location of the data
c
c       Local
c
        integer(hsize_t) :: h5dimsf(1:nd) ! file space dimensions
        integer(hsize_t) :: h5dimsm(1:nd) ! memory space dimensions
        integer(hsize_t) :: h5dimsv(1:nd) ! offset vector
        integer(hsize_t) :: h5dimss(1:nd) ! stride (file)
        integer(hsize_t) :: h5dimsc(1:nd) ! count: #blocks in each dim
        integer(hid_t)   :: h5fspaceid    ! file space ID
        integer(hid_t)   :: h5mspaceid    ! memory space ID
        integer(hid_t)   :: h5dsetid      ! data set ID
        integer(hid_t)   :: h5dmemtype    ! data memory type
        integer          :: h5select      ! selection parameter
c
        integer          :: ii            ! counter
c
        integer :: comm_size ! the size of communicator comm
        integer :: comm_rank ! the rank in communicator comm
        integer :: wrld_size ! the size of communicator MPI_COMM_WORLD
        integer :: wrld_rank ! the rank in communicator MPI_COMM_WORLD
        integer :: mpierr    ! the MPI error flag
c
        logical :: dataset_exists   ! does the data exist already?
        logical :: global_size_zero ! is the global dataset size 0?
        logical :: local_size_zero  ! is the local block size 0?
c
c       Code
c
        call mpi_comm_size(MPI_COMM_WORLD,wrld_size,mpierr)
        call mpi_comm_size(comm,comm_size,mpierr)
        call mpi_comm_rank(MPI_COMM_WORLD,wrld_rank,mpierr)
        call mpi_comm_rank(comm,comm_rank,mpierr)
c
        global_size_zero = .false.
        local_size_zero  = .false.
        do ii = 1, nd
          global_size_zero = global_size_zero.or.(ng(ii).le.0)
          local_size_zero  = local_size_zero .or.(nl(ii).le.0)
        enddo
        if (global_size_zero) then
          ! This dataset holds no data on any processor so just do
          ! nothing and return
          return
        endif
        h5dimsf = ng
        h5dimsm = nl
        h5dimsv = nv
        h5dimss = 1
        h5dimsc = 1
        h5dmemtype = H5T_NATIVE_DCMPLX
        h5select   = H5S_SELECT_SET_F
        call h5screate_simple_f(nd,h5dimsm,h5mspaceid,h5err)
        if (h5err.ne.0) then
          write(*,*)'ERROR: fread_zn: h5screate_simple_f failed',path
        endif
c
c       call h5lexists_f(fileid,path,dataset_exists,h5err)
        dataset_exists = fexist(fileid,path)
        if (.not.dataset_exists) then
          write(*,*)wrld_rank," dataset: ",path," does not exist!"
          call h5sclose_f(h5mspaceid,h5err)
          if (h5err.ne.0) then
            write(*,*)'ERROR: fread_zn: h5sclose_f failed',path
          endif
          return
        endif
        call h5dopen_f(fileid,path,h5dsetid,h5err)
        if (h5err.ne.0) then
          write(*,*)'ERROR: fread_zn: h5dopen_f failed',path
        endif
c
        call h5dget_space_f(h5dsetid,h5fspaceid,h5err)
        if (h5err.ne.0) then
          write(*,*)'ERROR: fread_zn: h5dget_space_f failed',path
        endif
        call h5sselect_hyperslab_f(h5fspaceid,h5select,h5dimsv,h5dimsc,
     &                             h5err,h5dimss,h5dimsm)
        if (h5err.ne.0) then
          write(*,*)'ERROR: fread_zn: h5sselect_hyperslab_f failed',path
        endif
c
        if (.not.local_size_zero) then
          call h5dread_f(h5dsetid,h5dmemtype,f_ptr,
     &                   h5err,mem_space_id=h5mspaceid,
     &                   file_space_id=h5fspaceid)
          if (h5err.ne.0) then
            write(*,*)'ERROR: fread_zn: h5dread_f failed',path
          endif
        else
          ! The local block size is zero, so no I/O can be done here.
          ! The other operations on the file check the consistency
          ! of the file across a number of processors. The data layout
          ! should be the same on all processors.
        endif
c
        call h5dclose_f(h5dsetid,h5err)
        if (h5err.ne.0) then
          write(*,*)'ERROR: fread_zn: h5dclose_f failed',path
        endif
        call h5sclose_f(h5fspaceid,h5err)
        if (h5err.ne.0) then
          write(*,*)'ERROR: fread_zn: h5sclose_f h5fspaceid failed',path
        endif
        call h5sclose_f(h5mspaceid,h5err)
        if (h5err.ne.0) then
          write(*,*)'ERROR: fread_zn: h5sclose_f h5mspaceid failed',path
        endif
c
        call mpi_barrier(comm,mpierr)
c
        end subroutine fread_zn
c
c       ----------------------------------------------------------------
c
c       Write an 0-dimensional double precision array
c
        subroutine fwrite_d0(comm,darray,path,fileid)
        use hdf5
        implicit none
        integer, intent(in) :: comm ! communicator for proc set
        integer(hid_t), intent(in) :: fileid ! the file identifier
c
        character(len=*), intent(in) :: path ! where to store
c
        double precision, intent(in) :: darray ! the data
c
c       Local
c
        double precision :: buf(1)
        integer :: nd    = 1
        integer :: ng(1) = 1
        integer :: nl(1) = 1
        integer :: nv(1) = 0
c
c       Code
c
        buf = darray
        call fwrite_dn(comm,buf,nd,ng,nl,nv,path,fileid)
c
        end subroutine fwrite_d0
c
c       ----------------------------------------------------------------
c
c       Write an 0-dimensional integer array
c
        subroutine fwrite_i0(comm,iarray,path,fileid)
        use hdf5
        implicit none
        integer, intent(in) :: comm ! communicator for proc set
        integer(hid_t), intent(in) :: fileid ! the file identifier
c
        character(len=*), intent(in) :: path ! where to store
c
        integer, intent(in) :: iarray ! the data
c
c       Local
c
        integer :: buf(1)
        integer :: nd    = 1
        integer :: ng(1) = 1
        integer :: nl(1) = 1
        integer :: nv(1) = 0
c
c       Code
c
        buf = iarray
        call fwrite_in(comm,buf,nd,ng,nl,nv,path,fileid)
c
        end subroutine fwrite_i0
c
c       ----------------------------------------------------------------
c
c       Write an 0-dimensional logical array
c
        subroutine fwrite_l0(comm,larray,path,fileid)
        use hdf5
        implicit none
        integer, intent(in) :: comm ! communicator for proc set
        integer(hid_t), intent(in) :: fileid ! the file identifier
c
        character(len=*), intent(in) :: path ! where to store
c
        logical, intent(in) :: larray ! the data
c
c       Local
c
        logical :: buf(1)
        integer :: nd    = 1
        integer :: ng(1) = 1
        integer :: nl(1) = 1
        integer :: nv(1) = 0
c
c       Code
c
        buf = larray
        call fwrite_ln(comm,buf,nd,ng,nl,nv,path,fileid)
c
        end subroutine fwrite_l0
c
c       ----------------------------------------------------------------
c
c       Write an 1-dimensional logical array
c
        subroutine fwrite_l1(comm,larray,ng,nv,path,fileid)
        use hdf5
        implicit none
        integer, intent(in) :: comm ! proc sub-set communicator
        integer(hid_t), intent(in) :: fileid ! the file identifier
c
        integer, intent(in) :: ng ! global dimension
        integer, intent(in) :: nv ! vector offset
c
        character(len=*), intent(in) :: path ! where to store
c
        logical, intent(in) :: larray(:) ! the data
c
c       Local
c
        integer :: nd = 1
        integer :: nnl(1)
        integer :: nng(1)
        integer :: nnv(1)
c
c       Code
c
        nnl(1) = size(larray,1)
        nng(1) = ng
        nnv(1) = nv
        call fwrite_ln(comm,larray,nd,nng,nnl,nnv,path,fileid)
c
        end subroutine fwrite_l1
c
c       ----------------------------------------------------------------
c
c       Write an 2-dimensional logical array
c
        subroutine fwrite_l2(comm,larray,ng,nv,path,fileid)
        use hdf5
        implicit none
        integer, intent(in) :: comm ! proc sub-set communicator
        integer(hid_t), intent(in) :: fileid ! the file identifier
c
        integer, intent(in) :: ng(2) ! global dimensions
        integer, intent(in) :: nv(2) ! vector offsets
c
        character(len=*), intent(in) :: path ! where to store
c
        logical, intent(in) :: larray(:,:) ! the data
c
c       Local
c
        integer :: nd = 2
        integer :: nl(2)
        integer :: nt
c
c       Code
c
        nl(1) = size(larray,1)
        nl(2) = size(larray,2)
        nt    = size(larray)
        call fwrite_ln(comm,reshape(larray,(/nt/)),nd,ng,nl,nv,path,
     &                 fileid)
c
        end subroutine fwrite_l2
c
c       ----------------------------------------------------------------
c
c       Write an 0-dimensional text array
c
        subroutine fwrite_t0(comm,tarray,path,fileid)
        use hdf5
        implicit none
        integer, intent(in) :: comm ! communicator for proc set
        integer(hid_t), intent(in) :: fileid ! the file identifier
c
        character(len=*), intent(in) :: path ! where to store
c
        character(len=*), intent(in), target :: tarray ! the data
c
c       Local
c
        integer(hsize_t) :: h5dims(1)
        integer(hsize_t) :: sizetxt
        integer(hid_t)   :: h5mtype  ! the memory type
        integer(hid_t)   :: h5ftype  ! the file type
        integer(hid_t)   :: h5dspaceid
        integer(hid_t)   :: h5dsetid
        integer(hid_t)   :: stringid
        type(c_ptr)      :: f_ptr
c
        integer :: mpierr
        integer :: irank
c
        logical :: dataset_exists ! does the data exist already?
c
c       Code
c
        call mpi_comm_rank(comm,irank,mpierr)
        h5dims = 1
        sizetxt = len(tarray)
        stringid = H5T_FORTRAN_S1
        call h5tcopy_f(stringid,h5mtype,h5err)
        call h5tcopy_f(stringid,h5ftype,h5err)
        call h5tset_size_f(h5mtype,sizetxt,h5err)
        call h5tset_size_f(h5ftype,sizetxt,h5err)
        dataset_exists = fexist(fileid,path)
        if (dataset_exists) then
          call h5dopen_f(fileid,path,h5dsetid,h5err)
        else
          call h5screate_simple_f(1,h5dims,h5dspaceid,h5err)
c         call h5dcreate_f(fileid,path,stringid,h5dspaceid,
c    &                     h5dsetid,h5err)
          call fcreate(fileid,path,h5ftype,h5dspaceid,h5dsetid)
          call h5sclose_f(h5dspaceid,h5err)
        endif
        if (irank.eq.0) then
          f_ptr = c_loc(tarray)
          call h5dwrite_f(h5dsetid,h5mtype,f_ptr,h5err)
        endif
        call h5dclose_f(h5dsetid,h5err)
        call h5tclose_f(h5mtype,h5err) 
        call h5tclose_f(h5ftype,h5err) 
c
        end subroutine fwrite_t0
c
c       ----------------------------------------------------------------
c
c       Write an 0-dimensional double complex array
c
        subroutine fwrite_z0(comm,zarray,path,fileid)
        use hdf5
        implicit none
        integer, intent(in) :: comm ! communicator for proc set
        integer(hid_t), intent(in) :: fileid ! the file identifier
c
        character(len=*), intent(in) :: path ! where to store
c
        double complex, intent(in), target :: zarray ! the data
c
c       Local
c
        double complex :: buf(1)
        integer :: nd    = 1
        integer :: ng(1) = 1
        integer :: nl(1) = 1
        integer :: nv(1) = 0
c
c       Code
c
        buf = zarray
        call fwrite_zn(comm,buf,nd,ng,nl,nv,path,fileid)
c
        end subroutine fwrite_z0
c
c       ----------------------------------------------------------------
c
c       Read an 0-dimensional double precision array
c
        subroutine fread_d0(comm,darray,path,fileid)
        use hdf5
        implicit none
        integer, intent(in) :: comm ! communicator for proc set
        integer(hid_t), intent(in) :: fileid ! the file identifier
c       
        character(len=*), intent(in) :: path ! where to store
c
        double precision, intent(out), target :: darray ! the data
c
c       Local
c
        type(c_ptr) :: f_ptr
        integer :: nd    = 1
        integer :: ng(1) = 1
        integer :: nl(1) = 1
        integer :: nv(1) = 0
c
c       Code
c
        f_ptr = c_loc(darray)
        call fread_dn(comm,f_ptr,nd,ng,nl,nv,path,fileid)
c
        end subroutine fread_d0
c
c       ----------------------------------------------------------------
c
c       Read an 0-dimensional integer array
c
        subroutine fread_i0(comm,iarray,path,fileid)
        use hdf5
        implicit none
        integer, intent(in) :: comm ! communicator for proc set
        integer(hid_t), intent(in) :: fileid ! the file identifier
c       
        character(len=*), intent(in) :: path ! where to store
c
        integer, intent(out), target :: iarray ! the data
c
c       Local
c
        integer :: nd    = 1
        integer :: ng(1) = 1
        integer :: nl(1) = 1
        integer :: nv(1) = 0
c
        type(c_ptr) :: f_ptr
c
c       Code
c
        f_ptr = c_loc(iarray)
        call fread_in(comm,f_ptr,nd,ng,nl,nv,path,fileid)
c
        end subroutine fread_i0
c
c       ----------------------------------------------------------------
c
c       Read an 0-dimensional logical array
c
        subroutine fread_l0(comm,larray,path,fileid)
        use hdf5
        implicit none
        integer, intent(in) :: comm ! communicator for proc set
        integer(hid_t), intent(in) :: fileid ! the file identifier
c       
        character(len=*), intent(in) :: path ! where to store
c
        logical, intent(out), target :: larray ! the data
c
c       Local
c
        integer :: nd    = 1
        integer :: ng(1) = 1
        integer :: nl(1) = 1
        integer :: nv(1) = 0
c
        type(c_ptr) :: f_ptr
c
c       Code
c
        f_ptr = c_loc(larray)
        call fread_ln(comm,f_ptr,nd,ng,nl,nv,path,fileid)
c
        end subroutine fread_l0
c
c       ----------------------------------------------------------------
c
c       Read an 1-dimensional logical array
c
        subroutine fread_l1(comm,larray,ng,nv,path,fileid)
        use hdf5
        implicit none
        integer, intent(in) :: comm ! proc sub-set communicator
        integer(hid_t), intent(in) :: fileid ! the file identifier
c
        integer, intent(in) :: ng ! global dimension
        integer, intent(in) :: nv ! vector offset
c       
        character(len=*), intent(in) :: path ! where to store
c
        logical, intent(out), target :: larray(:) ! the data
c
c       Local
c
        integer :: nd = 1
        integer :: nnl(1)
        integer :: nng(1)
        integer :: nnv(1)
c
        type(c_ptr) :: f_ptr
c
c       Code
c
        nnl(1)  = size(larray,1)
        nng(1)  = ng
        nnv(1)  = nv
        f_ptr = c_loc(larray)
        call fread_ln(comm,f_ptr,nd,nng,nnl,nnv,path,fileid)
c
        end subroutine fread_l1
c
c       ----------------------------------------------------------------
c
c       Read an 2-dimensional logical array
c
        subroutine fread_l2(comm,larray,ng,nv,path,fileid)
        use hdf5
        implicit none
        integer, intent(in) :: comm ! proc sub-set communicator
        integer(hid_t), intent(in) :: fileid ! the file identifier
c
        integer, intent(in) :: ng(2) ! global dimensions
        integer, intent(in) :: nv(2) ! vector offsets
c       
        character(len=*), intent(in) :: path ! where to store
c
        logical, intent(out), target :: larray(:,:) ! the data
c
c       Local
c
        integer :: nd = 2
        integer :: nl(2)
c
        type(c_ptr) :: f_ptr
c
c       Code
c
        nl(1) = size(larray,1)
        nl(2) = size(larray,2)
        f_ptr = c_loc(larray)
        call fread_dn(comm,f_ptr,nd,ng,nl,nv,path,fileid)
c
        end subroutine fread_l2
c
c       ----------------------------------------------------------------
c
c       Read an 0-dimensional text array
c
        subroutine fread_t0(comm,tarray,path,fileid)
        use hdf5
        implicit none
        integer, intent(in) :: comm ! communicator for proc set
        integer(hid_t), intent(in) :: fileid ! the file identifier
c
        character(len=*), intent(in) :: path ! where to store
c
        character(len=*), intent(out), target :: tarray ! the data
c
c       Local
c
        integer(hsize_t) :: h5dims(1)
        integer(hsize_t) :: sizetxt
        integer(hid_t)   :: h5mtype ! memory type
        integer(hid_t)   :: h5ftype ! file type
        integer(hid_t)   :: h5dspaceid
        integer(hid_t)   :: h5dsetid
        integer(hid_t)   :: stringid
        type(c_ptr)      :: f_ptr
c
        integer :: mpierr
        integer :: wrld_rank
c
        logical :: dataset_exists ! does the data exist already?
c
c       Code
c
        call mpi_comm_rank(MPI_COMM_WORLD,wrld_rank,mpierr)
        h5dims = 1
        sizetxt = len(tarray)
        stringid = H5T_FORTRAN_S1
        call h5tcopy_f(stringid,h5mtype,h5err)
        call h5tset_size_f(h5mtype,sizetxt,h5err)
c
c       call h5lexists_f(fileid,path,dataset_exists,h5err)
        dataset_exists = fexist(fileid,path)
        if (.not.dataset_exists) then
          write(*,*)wrld_rank," dataset: ",path," does not exist!"
c         call h5sclose_f(h5mspaceid,h5err)
          return
        endif
        call h5dopen_f(fileid,path,h5dsetid,h5err)
c       call h5dget_type_f(h5dsetid,h5ftype,h5err)
c       call h5dget_space_f(h5dsetid,h5dspaceid,h5err)
c
        f_ptr = c_loc(tarray)
        call h5dread_f(h5dsetid,h5mtype,f_ptr,h5err)
        call h5dclose_f(h5dsetid,h5err)
c       call h5sclose_f(h5dspaceid,h5err)
c       call h5tclose_f(h5ftype,h5err)
        call h5tclose_f(h5mtype,h5err)
c
        end subroutine fread_t0
c
c       ----------------------------------------------------------------
c
c       Read an 0-dimensional double complex array
c
        subroutine fread_z0(comm,zarray,path,fileid)
        use hdf5
        implicit none
        integer, intent(in) :: comm ! communicator for proc set
        integer(hid_t), intent(in) :: fileid ! the file identifier
c       
        character(len=*), intent(in) :: path ! where to store
c
        double complex, intent(out), target :: zarray ! the data
c
c       Local
c
        integer :: nd    = 1
        integer :: ng(1) = 1
        integer :: nl(1) = 1
        integer :: nv(1) = 0
c
        type(c_ptr) :: f_ptr
c
c       Code
c
        f_ptr = c_loc(zarray)
        call fread_zn(comm,f_ptr,nd,ng,nl,nv,path,fileid)
c
        end subroutine fread_z0
c
c       ----------------------------------------------------------------
c
c       Write an 1-dimensional double precision array
c
        subroutine fwrite_d1(comm,darray,ng,nv,path,fileid)
        use hdf5
        implicit none
        integer, intent(in) :: comm ! proc sub-set communicator
        integer(hid_t), intent(in) :: fileid ! the file identifier
c
        integer, intent(in) :: ng ! global dimension
        integer, intent(in) :: nv ! vector offset
c
        character(len=*), intent(in) :: path ! where to store
c
        double precision, intent(in) :: darray(:) ! the data
c
c       Local
c
        integer :: nd = 1
        integer :: nnl(1)
        integer :: nng(1)
        integer :: nnv(1)
c
c       Code
c
        nnl(1) = size(darray,1)
        nng(1) = ng
        nnv(1) = nv
        call fwrite_dn(comm,darray,nd,nng,nnl,nnv,path,fileid)
c
        end subroutine fwrite_d1
c
c       ----------------------------------------------------------------
c
c       Write an 1-dimensional integer array
c
        subroutine fwrite_i1(comm,iarray,ng,nv,path,fileid)
        use hdf5
        implicit none
        integer, intent(in) :: comm ! proc sub-set communicator
        integer(hid_t), intent(in) :: fileid ! the file identifier
c
        integer, intent(in) :: ng ! global dimension
        integer, intent(in) :: nv ! vector offset
c
        character(len=*), intent(in) :: path ! where to store
c
        integer, intent(in) :: iarray(:) ! the data
c
c       Local
c
        integer :: nd = 1
        integer :: nnl(1)
        integer :: nng(1)
        integer :: nnv(1)
c
c       Code
c
        nnl(1) = size(iarray,1)
        nng(1) = ng
        nnv(1) = nv
        call fwrite_in(comm,iarray,nd,nng,nnl,nnv,path,fileid)
c
        end subroutine fwrite_i1
c
c       ----------------------------------------------------------------
c
c       Write an 1-dimensional character array
c
        subroutine fwrite_t1(comm,tarray,ng,nv,path,fileid)
        use hdf5
        implicit none
        integer, intent(in) :: comm ! proc sub-set communicator
        integer(hid_t), intent(in) :: fileid ! the file identifier
c
        integer, intent(in) :: ng ! global dimension
        integer, intent(in) :: nv ! vector offset
c
        character(len=*), intent(in) :: path ! where to store
c
        character(len=*), intent(in) :: tarray(:) ! the data
c
c       Local
c
        integer :: nd = 1
        integer :: nnl(1)
        integer :: nng(1)
        integer :: nnv(1)
c
c       Code
c
        nnl(1) = size(tarray,1)
        nng(1) = ng
        nnv(1) = nv
        call fwrite_tn(comm,tarray,nd,nng,nnl,nnv,path,fileid)
c
        end subroutine fwrite_t1
c
c       ----------------------------------------------------------------
c
c       Write an 1-dimensional double complex array
c
        subroutine fwrite_z1(comm,zarray,ng,nv,path,fileid)
        use hdf5
        implicit none
        integer, intent(in) :: comm ! proc sub-set communicator
        integer(hid_t), intent(in) :: fileid ! the file identifier
c
        integer, intent(in) :: ng ! global dimension
        integer, intent(in) :: nv ! vector offset
c
        character(len=*), intent(in) :: path ! where to store
c
        double complex, intent(in), target :: zarray(:) ! the data
c
c       Local
c
        integer :: nd = 1
        integer :: nnl(1)
        integer :: nng(1)
        integer :: nnv(1)
c
c       Code
c
        nnl(1) = size(zarray,1)
        nng(1) = ng
        nnv(1) = nv
        call fwrite_zn(comm,zarray,nd,nng,nnl,nnv,path,fileid)
c
        end subroutine fwrite_z1
c
c       ----------------------------------------------------------------
c
c       Read an 1-dimensional double precision array
c
        subroutine fread_d1(comm,darray,ng,nv,path,fileid)
        use hdf5
        implicit none
        integer, intent(in) :: comm ! proc sub-set communicator
        integer(hid_t), intent(in) :: fileid ! the file identifier
c
        integer, intent(in) :: ng ! global dimension
        integer, intent(in) :: nv ! vector offset
c       
        character(len=*), intent(in) :: path ! where to store
c
        double precision, intent(out), target :: darray(:) ! the data
c
c       Local
c
        integer :: nd = 1
        integer :: nnl(1)
        integer :: nng(1)
        integer :: nnv(1)
c
        type(c_ptr) :: f_ptr
c
c       Code
c
        nnl(1)  = size(darray,1)
        nng(1)  = ng
        nnv(1)  = nv
        f_ptr = c_loc(darray)
        call fread_dn(comm,f_ptr,nd,nng,nnl,nnv,path,fileid)
c
        end subroutine fread_d1
c
c       ----------------------------------------------------------------
c
c       Read an 1-dimensional integer array
c
        subroutine fread_i1(comm,iarray,ng,nv,path,fileid)
        use hdf5
        implicit none
        integer, intent(in) :: comm ! proc sub-set communicator
        integer(hid_t), intent(in) :: fileid ! the file identifier
c
        integer, intent(in) :: ng ! global dimension
        integer, intent(in) :: nv ! vector offset
c       
        character(len=*), intent(in) :: path ! where to store
c
        integer, intent(out), target :: iarray(:) ! the data
c
c       Local
c
        integer :: nd = 1
        integer :: nnl(1)
        integer :: nng(1)
        integer :: nnv(1)
c
        type(c_ptr) :: f_ptr
c
c       Code
c
        nnl(1)  = size(iarray,1)
        nng(1)  = ng
        nnv(1)  = nv
        f_ptr = c_loc(iarray)
        call fread_in(comm,f_ptr,nd,nng,nnl,nnv,path,fileid)
c
        end subroutine fread_i1
c
c       ----------------------------------------------------------------
c
c       Read an 1-dimensional character array
c
        subroutine fread_t1(comm,tarray,ng,nv,path,fileid)
        use hdf5
        implicit none
        integer, intent(in) :: comm ! proc sub-set communicator
        integer(hid_t), intent(in) :: fileid ! the file identifier
c
        integer, intent(in) :: ng ! global dimension
        integer, intent(in) :: nv ! vector offset
c       
        character(len=*), intent(in) :: path ! where to store
c
        character(len=*), intent(out), target :: tarray(:) ! the data
c
c       Local
c
        integer :: nd = 1
        integer :: nnl(1)
        integer :: nng(1)
        integer :: nnv(1)
        integer :: nnt
c
        type(c_ptr) :: f_ptr
c
c       Code
c
        nnl(1)  = size(tarray,1)
        nng(1)  = ng
        nnv(1)  = nv
        nnt     = len(tarray(1))
        f_ptr = c_loc(tarray)
        call fread_tn(comm,f_ptr,nnt,nd,nng,nnl,nnv,path,fileid)
c
        end subroutine fread_t1
c
c       ----------------------------------------------------------------
c
c       Read an 1-dimensional double complex array
c
        subroutine fread_z1(comm,zarray,ng,nv,path,fileid)
        use hdf5
        implicit none
        integer, intent(in) :: comm ! proc sub-set communicator
        integer(hid_t), intent(in) :: fileid ! the file identifier
c
        integer, intent(in) :: ng ! global dimension
        integer, intent(in) :: nv ! vector offset
c       
        character(len=*), intent(in) :: path ! where to store
c
        double complex, intent(out), target :: zarray(:) ! the data
c
c       Local
c
        integer :: nd = 1
        integer :: nnl(1)
        integer :: nng(1)
        integer :: nnv(1)
c
        type(c_ptr) :: f_ptr
c
c       Code
c
        nnl(1) = size(zarray,1)
        nng(1) = ng
        nnv(1) = nv
        f_ptr  = c_loc(zarray)
        call fread_zn(comm,f_ptr,nd,nng,nnl,nnv,path,fileid)
c
        end subroutine fread_z1
c
c       ----------------------------------------------------------------
c
c       Write an 2-dimensional double precision array
c
        subroutine fwrite_d2(comm,darray,ng,nv,path,fileid)
        use hdf5
        implicit none
        integer, intent(in) :: comm ! proc sub-set communicator
        integer(hid_t), intent(in) :: fileid ! the file identifier
c
        integer, intent(in) :: ng(2) ! global dimensions
        integer, intent(in) :: nv(2) ! vector offsets
c
        character(len=*), intent(in) :: path ! where to store
c
        double precision, intent(in) :: darray(:,:) ! the data
c
c       Local
c
        integer :: nd = 2
        integer :: nl(2)
        integer :: nt
c
c       Code
c
        nl(1) = size(darray,1)
        nl(2) = size(darray,2)
        nt    = size(darray)
        call fwrite_dn(comm,reshape(darray,(/nt/)),nd,ng,nl,nv,path,
     &                 fileid)
c
        end subroutine fwrite_d2
c
c       ----------------------------------------------------------------
c
c       Write an 2-dimensional integer array
c
        subroutine fwrite_i2(comm,iarray,ng,nv,path,fileid)
        use hdf5
        implicit none
        integer, intent(in) :: comm ! proc sub-set communicator
        integer(hid_t), intent(in) :: fileid ! the file identifier
c
        integer, intent(in) :: ng(2) ! global dimensions
        integer, intent(in) :: nv(2) ! vector offsets
c
        character(len=*), intent(in) :: path ! where to store
c
        integer, intent(in) :: iarray(:,:) ! the data
c
c       Local
c
        integer :: nd = 2
        integer :: nl(2)
        integer :: nt
c
c       Code
c
        nl(1) = size(iarray,1)
        nl(2) = size(iarray,2)
        nt    = size(iarray)
        call fwrite_in(comm,reshape(iarray,(/nt/)),nd,ng,nl,nv,path,
     &                 fileid)
c
        end subroutine fwrite_i2
c
c       ----------------------------------------------------------------
c
c       Write an 2-dimensional character array
c
        subroutine fwrite_t2(comm,tarray,ng,nv,path,fileid)
        use hdf5
        implicit none
        integer, intent(in) :: comm ! proc sub-set communicator
        integer(hid_t), intent(in) :: fileid ! the file identifier
c
        integer, intent(in) :: ng(2) ! global dimension
        integer, intent(in) :: nv(2) ! vector offset
c
        character(len=*), intent(in) :: path ! where to store
c
        character(len=*), intent(in) :: tarray(:,:) ! the data
c
c       Local
c
        integer :: nd = 2
        integer :: nl(2)
        integer :: nt
c
c       Code
c
        nl(1) = size(tarray,1)
        nl(2) = size(tarray,2)
        nt    = size(tarray)
        call fwrite_tn(comm,reshape(tarray,(/nt/)),nd,ng,nl,nv,path,
     &                 fileid)
c
        end subroutine fwrite_t2
c
c       ----------------------------------------------------------------
c
c       Write an 2-dimensional double complex array
c
        subroutine fwrite_z2(comm,zarray,ng,nv,path,fileid)
        use hdf5
        implicit none
        integer, intent(in) :: comm ! proc sub-set communicator
        integer(hid_t), intent(in) :: fileid ! the file identifier
c
        integer, intent(in) :: ng(2) ! global dimensions
        integer, intent(in) :: nv(2) ! vector offsets
c
        character(len=*), intent(in) :: path ! where to store
c
        double complex, intent(in), target :: zarray(:,:) ! the data
c
c       Local
c
        integer :: nd = 2
        integer :: nl(2)
        integer :: nt
c
c       Code
c
        nl(1) = size(zarray,1)
        nl(2) = size(zarray,2)
        nt    = size(zarray)
        call fwrite_zn(comm,reshape(zarray,(/nt/)),nd,ng,nl,nv,path,
     &                 fileid)
c
        end subroutine fwrite_z2
c
c       ----------------------------------------------------------------
c
c       Read an 2-dimensional double precision array
c
        subroutine fread_d2(comm,darray,ng,nv,path,fileid)
        use hdf5
        implicit none
        integer, intent(in) :: comm ! proc sub-set communicator
        integer(hid_t), intent(in) :: fileid ! the file identifier
c
        integer, intent(in) :: ng(2) ! global dimensions
        integer, intent(in) :: nv(2) ! vector offsets
c       
        character(len=*), intent(in) :: path ! where to store
c
        double precision, intent(out), target :: darray(:,:) ! the data
c
c       Local
c
        integer :: nd = 2
        integer :: nl(2)
c
        type(c_ptr) :: f_ptr
c
c       Code
c
        nl(1) = size(darray,1)
        nl(2) = size(darray,2)
        f_ptr = c_loc(darray)
        call fread_dn(comm,f_ptr,nd,ng,nl,nv,path,fileid)
c
        end subroutine fread_d2
c
c       ----------------------------------------------------------------
c
c       Read an 2-dimensional integer array
c
        subroutine fread_i2(comm,iarray,ng,nv,path,fileid)
        use hdf5
        implicit none
        integer, intent(in) :: comm ! proc sub-set communicator
        integer(hid_t), intent(in) :: fileid ! the file identifier
c
        integer, intent(in) :: ng(2) ! global dimensions
        integer, intent(in) :: nv(2) ! vector offsets
c
        character(len=*), intent(in) :: path ! where to store
c
        integer, intent(out), target :: iarray(:,:) ! the data
c
c       Local
c
        integer :: nd = 2
        integer :: nl(2)
c
        type(c_ptr) :: f_ptr
c
c       Code
c
        nl(1) = size(iarray,1)
        nl(2) = size(iarray,2)
        f_ptr = c_loc(iarray)
        call fread_in(comm,f_ptr,nd,ng,nl,nv,path,fileid)
c
        end subroutine fread_i2
c
c       ----------------------------------------------------------------
c
c       Read an 2-dimensional character array
c
        subroutine fread_t2(comm,tarray,ng,nv,path,fileid)
        use hdf5
        implicit none
        integer, intent(in) :: comm ! proc sub-set communicator
        integer(hid_t), intent(in) :: fileid ! the file identifier
c
        integer, intent(in) :: ng(2) ! global dimension
        integer, intent(in) :: nv(2) ! vector offset
c       
        character(len=*), intent(in) :: path ! where to store
c
        character(len=*), intent(out), target :: tarray(:,:) ! the data
c
c       Local
c
        integer :: nd = 2
        integer :: nl(2)
        integer :: nt
c
        type(c_ptr) :: f_ptr
c
c       Code
c
        nl(1) = size(tarray,1)
        nl(2) = size(tarray,2)
        nt    = len(tarray(1,1))
        f_ptr = c_loc(tarray)
        call fread_tn(comm,f_ptr,nt,nd,ng,nl,nv,path,fileid)
c
        end subroutine fread_t2
c
c       ----------------------------------------------------------------
c
c       Read an 2-dimensional double complex array
c
        subroutine fread_z2(comm,zarray,ng,nv,path,fileid)
        use hdf5
        implicit none
        integer, intent(in) :: comm ! proc sub-set communicator
        integer(hid_t), intent(in) :: fileid ! the file identifier
c
        integer, intent(in) :: ng(2) ! global dimensions
        integer, intent(in) :: nv(2) ! vector offsets
c       
        character(len=*), intent(in) :: path ! where to store
c
        double complex, intent(out), target :: zarray(:,:) ! the data
c
c       Local
c
        integer :: nd = 2
        integer :: nl(2)
c
        type(c_ptr) :: f_ptr
c
c       Code
c
        nl(1) = size(zarray,1)
        nl(2) = size(zarray,2)
        f_ptr = c_loc(zarray)
        call fread_zn(comm,f_ptr,nd,ng,nl,nv,path,fileid)
c
        end subroutine fread_z2
c
c       ----------------------------------------------------------------
c
c       Write an 3-dimensional double precision array
c
        subroutine fwrite_d3(comm,darray,ng,nv,path,fileid)
        use hdf5
        implicit none
        integer, intent(in) :: comm ! proc sub-set communicator
        integer(hid_t), intent(in) :: fileid ! the file identifier
c
        integer, intent(in) :: ng(3) ! global dimensions
        integer, intent(in) :: nv(3) ! vector offsets
c
        character(len=*), intent(in) :: path ! where to store
c
        double precision, intent(in) :: darray(:,:,:) ! the data
c
c       Local
c
        integer :: nd = 3
        integer :: nl(3)
        integer :: nt
c
c       Code
c
        nl(1) = size(darray,1)
        nl(2) = size(darray,2)
        nl(3) = size(darray,3)
        nt    = size(darray)
        call fwrite_dn(comm,reshape(darray,(/nt/)),nd,ng,nl,nv,path,
     &                 fileid)
c
        end subroutine fwrite_d3
c
c       ----------------------------------------------------------------
c
c       Write an 3-dimensional integer array
c
        subroutine fwrite_i3(comm,iarray,ng,nv,path,fileid)
        use hdf5
        implicit none
        integer, intent(in) :: comm ! proc sub-set communicator
        integer(hid_t), intent(in) :: fileid ! the file identifier
c
        integer, intent(in) :: ng(3) ! global dimensions
        integer, intent(in) :: nv(3) ! vector offsets
c
        character(len=*), intent(in) :: path ! where to store
c
        integer, intent(in) :: iarray(:,:,:) ! the data
c
c       Local
c
        integer :: nd = 3
        integer :: nl(3)
        integer :: nt
c
c       Code
c
        nl(1) = size(iarray,1)
        nl(2) = size(iarray,2)
        nl(3) = size(iarray,3)
        nt    = size(iarray)
        call fwrite_in(comm,reshape(iarray,(/nt/)),nd,ng,nl,nv,path,
     &                 fileid)
c
        end subroutine fwrite_i3
c
c       ----------------------------------------------------------------
c
c       Write an 3-dimensional character array
c
        subroutine fwrite_t3(comm,tarray,ng,nv,path,fileid)
        use hdf5
        implicit none
        integer, intent(in) :: comm ! proc sub-set communicator
        integer(hid_t), intent(in) :: fileid ! the file identifier
c
        integer, intent(in) :: ng(3) ! global dimension
        integer, intent(in) :: nv(3) ! vector offset
c
        character(len=*), intent(in) :: path ! where to store
c
        character(len=*), intent(in) :: tarray(:,:,:) ! the data
c
c       Local
c
        integer :: nd = 3
        integer :: nl(3)
        integer :: nt
c
c       Code
c
        nl(1) = size(tarray,1)
        nl(2) = size(tarray,2)
        nl(3) = size(tarray,3)
        nt    = size(tarray)
        call fwrite_tn(comm,reshape(tarray,(/nt/)),nd,ng,nl,nv,path,
     &                 fileid)
c
        end subroutine fwrite_t3
c
c       ----------------------------------------------------------------
c
c       Write an 3-dimensional double complex array
c
        subroutine fwrite_z3(comm,zarray,ng,nv,path,fileid)
        use hdf5
        implicit none
        integer, intent(in) :: comm ! proc sub-set communicator
        integer(hid_t), intent(in) :: fileid ! the file identifier
c
        integer, intent(in) :: ng(3) ! global dimensions
        integer, intent(in) :: nv(3) ! vector offsets
c
        character(len=*), intent(in) :: path ! where to store
c
        double complex, intent(in), target :: zarray(:,:,:) ! the data
c
c       Local
c
        integer :: nd = 3
        integer :: nl(3)
        integer :: nt
c
c       Code
c
        nl(1) = size(zarray,1)
        nl(2) = size(zarray,2)
        nl(3) = size(zarray,3)
        nt    = size(zarray)
        call fwrite_zn(comm,reshape(zarray,(/nt/)),nd,ng,nl,nv,path,
     &                 fileid)
c
        end subroutine fwrite_z3
c
c       ----------------------------------------------------------------
c
c       Read an 3-dimensional double precision array
c
        subroutine fread_d3(comm,darray,ng,nv,path,fileid)
        use hdf5
        implicit none
        integer, intent(in) :: comm ! proc sub-set communicator
        integer(hid_t), intent(in) :: fileid ! the file identifier
c
        integer, intent(in) :: ng(3) ! global dimensions
        integer, intent(in) :: nv(3) ! vector offsets
c
        character(len=*), intent(in) :: path ! where to store
c
        double precision, intent(out), target :: darray(:,:,:) ! the data
c
c       Local
c
        integer :: nd = 3
        integer :: nl(3)
c
        type(c_ptr) :: f_ptr
c
c       Code
c
        nl(1) = size(darray,1)
        nl(2) = size(darray,2)
        nl(3) = size(darray,3)
        f_ptr = c_loc(darray)
        call fread_dn(comm,f_ptr,nd,ng,nl,nv,path,fileid)
c
        end subroutine fread_d3
c
c       ----------------------------------------------------------------
c
c       Read an 3-dimensional integer array
c
        subroutine fread_i3(comm,iarray,ng,nv,path,fileid)
        use hdf5
        implicit none
        integer, intent(in) :: comm ! proc sub-set communicator
        integer(hid_t), intent(in) :: fileid ! the file identifier
c
        integer, intent(in) :: ng(3) ! global dimensions
        integer, intent(in) :: nv(3) ! vector offsets
c       
        character(len=*), intent(in) :: path ! where to store
c
        integer, intent(out), target :: iarray(:,:,:) ! the data
c
c       Local
c
        integer :: nd = 3
        integer :: nl(3)
c
        type(c_ptr) :: f_ptr
c
c       Code
c
        nl(1) = size(iarray,1)
        nl(2) = size(iarray,2)
        nl(3) = size(iarray,3)
        f_ptr = c_loc(iarray)
        call fread_in(comm,f_ptr,nd,ng,nl,nv,path,fileid)
c
        end subroutine fread_i3
c
c       ----------------------------------------------------------------
c
c       Read an 3-dimensional character array
c
        subroutine fread_t3(comm,tarray,ng,nv,path,fileid)
        use hdf5
        implicit none
        integer, intent(in) :: comm ! proc sub-set communicator
        integer(hid_t), intent(in) :: fileid ! the file identifier
c
        integer, intent(in) :: ng(3) ! global dimension
        integer, intent(in) :: nv(3) ! vector offset
c       
        character(len=*), intent(in) :: path ! where to store
c
        character(len=*), intent(out), target :: tarray(:,:,:) ! data
c
c       Local
c
        integer :: nd = 3
        integer :: nl(3)
        integer :: nt
c
        type(c_ptr) :: f_ptr
c
c       Code
c
        nl(1) = size(tarray,1)
        nl(2) = size(tarray,2)
        nl(3) = size(tarray,3)
        nt    = len(tarray(1,1,1))
        f_ptr = c_loc(tarray)
        call fread_tn(comm,f_ptr,nt,nd,ng,nl,nv,path,fileid)
c
        end subroutine fread_t3
c
c       ----------------------------------------------------------------
c
c       Read an 3-dimensional double complex array
c
        subroutine fread_z3(comm,zarray,ng,nv,path,fileid)
        use hdf5
        implicit none
        integer, intent(in) :: comm ! proc sub-set communicator
        integer(hid_t), intent(in) :: fileid ! the file identifier
c
        integer, intent(in) :: ng(3) ! global dimensions
        integer, intent(in) :: nv(3) ! vector offsets
c       
        character(len=*), intent(in) :: path ! where to store
c
        double complex, intent(out), target :: zarray(:,:,:) ! the data
c
c       Local
c
        integer :: nd = 3
        integer :: nl(3)
c
        type(c_ptr) :: f_ptr
c
c       Code
c
        nl(1) = size(zarray,1)
        nl(2) = size(zarray,2)
        nl(3) = size(zarray,3)
        f_ptr = c_loc(zarray)
        call fread_zn(comm,f_ptr,nd,ng,nl,nv,path,fileid)
c
        end subroutine fread_z3
c
c       ----------------------------------------------------------------
c
c       Write an 4-dimensional double precision array
c
        subroutine fwrite_d4(comm,darray,ng,nv,path,fileid)
        use hdf5
        implicit none
        integer, intent(in) :: comm ! proc sub-set communicator
        integer(hid_t), intent(in) :: fileid ! the file identifier
c
        integer, intent(in) :: ng(4) ! global dimensions
        integer, intent(in) :: nv(4) ! vector offsets
c
        character(len=*), intent(in) :: path ! where to store
c
        double precision, intent(in) :: darray(:,:,:,:) ! the data
c
c       Local
c
        integer :: nd = 4
        integer :: nl(4)
        integer :: nt
c
c       Code
c
        nl(1) = size(darray,1)
        nl(2) = size(darray,2)
        nl(3) = size(darray,3)
        nl(4) = size(darray,4)
        nt    = size(darray)
        call fwrite_dn(comm,reshape(darray,(/nt/)),nd,ng,nl,nv,path,
     &                 fileid)
c
        end subroutine fwrite_d4
c
c       ----------------------------------------------------------------
c
c       Write a 5-dimensional double precision array
c
        subroutine fwrite_d5(comm,darray,ng,nv,path,fileid)
        use hdf5
        implicit none
        integer, intent(in) :: comm ! proc sub-set communicator
        integer(hid_t), intent(in) :: fileid ! the file identifier
c
        integer, intent(in) :: ng(5) ! global dimensions
        integer, intent(in) :: nv(5) ! vector offsets
c
        character(len=*), intent(in) :: path ! where to store
c
        double precision, intent(in) :: darray(:,:,:,:,:) ! the data
c
c       Local
c
        integer :: nd = 5
        integer :: nl(5)
        integer :: nt
c
c       Code
c
        nl(1) = size(darray,1)
        nl(2) = size(darray,2)
        nl(3) = size(darray,3)
        nl(4) = size(darray,4)
        nl(5) = size(darray,5)
        nt    = size(darray)
        call fwrite_dn(comm,reshape(darray,(/nt/)),nd,ng,nl,nv,path,
     &                 fileid)
c
        end subroutine fwrite_d5
c
c       ----------------------------------------------------------------
c
c       Write a 6-dimensional double precision array
c
        subroutine fwrite_d6(comm,darray,ng,nv,path,fileid)
        use hdf5
        implicit none
        integer, intent(in) :: comm ! proc sub-set communicator
        integer(hid_t), intent(in) :: fileid ! the file identifier
c
        integer, intent(in) :: ng(6) ! global dimensions
        integer, intent(in) :: nv(6) ! vector offsets
c
        character(len=*), intent(in) :: path ! where to store
c
        double precision, intent(in) :: darray(:,:,:,:,:,:) ! the data
c
c       Local
c
        integer :: nd = 6
        integer :: nl(6)
        integer :: nt
c
c       Code
c
        nl(1) = size(darray,1)
        nl(2) = size(darray,2)
        nl(3) = size(darray,3)
        nl(4) = size(darray,4)
        nl(5) = size(darray,5)
        nl(6) = size(darray,6)
        nt    = size(darray)
        call fwrite_dn(comm,reshape(darray,(/nt/)),nd,ng,nl,nv,path,
     &                 fileid)
c
        end subroutine fwrite_d6
c
c       ----------------------------------------------------------------
c
c       Write an 7-dimensional double precision array
c
        subroutine fwrite_d7(comm,darray,ng,nv,path,fileid)
        use hdf5
        implicit none
        integer, intent(in) :: comm ! proc sub-set communicator
        integer(hid_t), intent(in) :: fileid ! the file identifier
c
        integer, intent(in) :: ng(7) ! global dimensions
        integer, intent(in) :: nv(7) ! vector offsets
c
        character(len=*), intent(in) :: path ! where to store
c
        double precision, intent(in) :: darray(:,:,:,:,:,:,:) ! the data
c
c       Local
c
        integer :: nd = 7
        integer :: nl(7)
        integer :: nt
c
c       Code
c
        nl(1) = size(darray,1)
        nl(2) = size(darray,2)
        nl(3) = size(darray,3)
        nl(4) = size(darray,4)
        nl(5) = size(darray,5)
        nl(6) = size(darray,6)
        nl(7) = size(darray,7)
        nt    = size(darray)
        call fwrite_dn(comm,reshape(darray,(/nt/)),nd,ng,nl,nv,path,
     &                 fileid)
c
        end subroutine fwrite_d7
c
c       ----------------------------------------------------------------
c
c       Write an 4-dimensional integer array
c
        subroutine fwrite_i4(comm,iarray,ng,nv,path,fileid)
        use hdf5
        implicit none
        integer, intent(in) :: comm ! proc sub-set communicator
        integer(hid_t), intent(in) :: fileid ! the file identifier
c
        integer, intent(in) :: ng(4) ! global dimensions
        integer, intent(in) :: nv(4) ! vector offsets
c
        character(len=*), intent(in) :: path ! where to store
c
        integer, intent(in) :: iarray(:,:,:,:) ! the data
c
c       Local
c
        integer :: nd = 4
        integer :: nl(4)
        integer :: nt
c
c       Code
c
        nl(1) = size(iarray,1)
        nl(2) = size(iarray,2)
        nl(3) = size(iarray,3)
        nl(4) = size(iarray,4)
        nt    = size(iarray)
        call fwrite_in(comm,reshape(iarray,(/nt/)),nd,ng,nl,nv,path,
     &                 fileid)
c
        end subroutine fwrite_i4
c
c       ----------------------------------------------------------------
c
c       Write an 4-dimensional double complex array
c
        subroutine fwrite_z4(comm,zarray,ng,nv,path,fileid)
        use hdf5
        implicit none
        integer, intent(in) :: comm ! proc sub-set communicator
        integer(hid_t), intent(in) :: fileid ! the file identifier
c
        integer, intent(in) :: ng(4) ! global dimensions
        integer, intent(in) :: nv(4) ! vector offsets
c
        character(len=*), intent(in) :: path ! where to store
c
        double complex, intent(in), target :: zarray(:,:,:,:) ! the data
c
c       Local
c
        integer :: nd = 4
        integer :: nl(4)
        integer :: nt
c
c       Code
c
        nl(1) = size(zarray,1)
        nl(2) = size(zarray,2)
        nl(3) = size(zarray,3)
        nl(4) = size(zarray,4)
        nt    = size(zarray)
        call fwrite_zn(comm,reshape(zarray,(/nt/)),nd,ng,nl,nv,path,
     &                 fileid)
c
        end subroutine fwrite_z4
c
c       ----------------------------------------------------------------
c
c       Read an 4-dimensional double precision array
c
        subroutine fread_d4(comm,darray,ng,nv,path,fileid)
        use hdf5
        implicit none
        integer, intent(in) :: comm ! proc sub-set communicator
        integer(hid_t), intent(in) :: fileid ! the file identifier
c
        integer, intent(in) :: ng(4) ! global dimensions
        integer, intent(in) :: nv(4) ! vector offsets
c       
        character(len=*), intent(in) :: path ! where to store
c
        double precision, intent(out), target :: darray(:,:,:,:) ! the data
c
c       Local
c
        integer :: nd = 4
        integer :: nl(4)
c
        type(c_ptr) :: f_ptr
c
c       Code
c
        nl(1) = size(darray,1)
        nl(2) = size(darray,2)
        nl(3) = size(darray,3)
        nl(4) = size(darray,4)
        f_ptr = c_loc(darray)
        call fread_dn(comm,f_ptr,nd,ng,nl,nv,path,fileid)
c
        end subroutine fread_d4
c
c       ----------------------------------------------------------------
c
c       Read a 5-dimensional double precision array
c
        subroutine fread_d5(comm,darray,ng,nv,path,fileid)
        use hdf5
        implicit none
        integer, intent(in) :: comm ! proc sub-set communicator
        integer(hid_t), intent(in) :: fileid ! the file identifier
c
        integer, intent(in) :: ng(5) ! global dimensions
        integer, intent(in) :: nv(5) ! vector offsets
c       
        character(len=*), intent(in) :: path ! where to store
c
        double precision, intent(out), target :: darray(:,:,:,:,:) ! the data
c
c       Local
c
        integer :: nd = 5
        integer :: nl(5)
c
        type(c_ptr) :: f_ptr
c
c       Code
c
        nl(1) = size(darray,1)
        nl(2) = size(darray,2)
        nl(3) = size(darray,3)
        nl(4) = size(darray,4)
        nl(5) = size(darray,5)
        f_ptr = c_loc(darray)
        call fread_dn(comm,f_ptr,nd,ng,nl,nv,path,fileid)
c
        end subroutine fread_d5
c
c       ----------------------------------------------------------------
c
c       Read a 6-dimensional double precision array
c
        subroutine fread_d6(comm,darray,ng,nv,path,fileid)
        use hdf5
        implicit none
        integer, intent(in) :: comm ! proc sub-set communicator
        integer(hid_t), intent(in) :: fileid ! the file identifier
c
        integer, intent(in) :: ng(6) ! global dimensions
        integer, intent(in) :: nv(6) ! vector offsets
c       
        character(len=*), intent(in) :: path ! where to store
c
        double precision, intent(out), target :: darray(:,:,:,:,:,:) ! the data
c
c       Local
c
        integer :: nd = 6
        integer :: nl(6)
c
        type(c_ptr) :: f_ptr
c
c       Code
c
        nl(1) = size(darray,1)
        nl(2) = size(darray,2)
        nl(3) = size(darray,3)
        nl(4) = size(darray,4)
        nl(5) = size(darray,5)
        nl(6) = size(darray,6)
        f_ptr = c_loc(darray)
        call fread_dn(comm,f_ptr,nd,ng,nl,nv,path,fileid)
c
        end subroutine fread_d6
c
c       ----------------------------------------------------------------
c
c       Read a 7-dimensional double precision array
c
        subroutine fread_d7(comm,darray,ng,nv,path,fileid)
        use hdf5
        implicit none
        integer, intent(in) :: comm ! proc sub-set communicator
        integer(hid_t), intent(in) :: fileid ! the file identifier
c
        integer, intent(in) :: ng(7) ! global dimensions
        integer, intent(in) :: nv(7) ! vector offsets
c       
        character(len=*), intent(in) :: path ! where to store
c
        double precision, intent(out), target :: darray(:,:,:,:,:,:,:) ! the data
c
c       Local
c
        integer :: nd = 7
        integer :: nl(7)
c
        type(c_ptr) :: f_ptr
c
c       Code
c
        nl(1) = size(darray,1)
        nl(2) = size(darray,2)
        nl(3) = size(darray,3)
        nl(4) = size(darray,4)
        nl(5) = size(darray,5)
        nl(6) = size(darray,6)
        nl(7) = size(darray,7)
        f_ptr = c_loc(darray)
        call fread_dn(comm,f_ptr,nd,ng,nl,nv,path,fileid)
c
        end subroutine fread_d7
c
c       ----------------------------------------------------------------
c
c       Read an 4-dimensional integer array
c
        subroutine fread_i4(comm,iarray,ng,nv,path,fileid)
        use hdf5
        implicit none
        integer, intent(in) :: comm ! proc sub-set communicator
        integer(hid_t), intent(in) :: fileid ! the file identifier
c
        integer, intent(in) :: ng(4) ! global dimensions
        integer, intent(in) :: nv(4) ! vector offsets
c
        character(len=*), intent(in) :: path ! where to store
c
        integer, intent(out), target :: iarray(:,:,:,:) ! the data
c
c       Local
c
        integer :: nd = 4
        integer :: nl(4)
c
        type(c_ptr) :: f_ptr
c
c       Code
c
        nl(1)  = size(iarray,1)
        nl(2)  = size(iarray,2)
        nl(3)  = size(iarray,3)
        nl(4)  = size(iarray,4)
        f_ptr = c_loc(iarray)
        call fread_in(comm,f_ptr,nd,ng,nl,nv,path,fileid)
c
        end subroutine fread_i4
c
c       ----------------------------------------------------------------
c
c       Read an 4-dimensional double complex array
c
        subroutine fread_z4(comm,zarray,ng,nv,path,fileid)
        use hdf5
        implicit none
        integer, intent(in) :: comm ! proc sub-set communicator
        integer(hid_t), intent(in) :: fileid ! the file identifier
c
        integer, intent(in) :: ng(4) ! global dimensions
        integer, intent(in) :: nv(4) ! vector offsets
c       
        character(len=*), intent(in) :: path ! where to store
c
        double complex, intent(out), target :: zarray(:,:,:,:) ! the data
c
c       Local
c
        integer :: nd = 4
        integer :: nl(4)
c
        type(c_ptr) :: f_ptr
c
c       Code
c
        nl(1) = size(zarray,1)
        nl(2) = size(zarray,2)
        nl(3) = size(zarray,3)
        nl(4) = size(zarray,4)
        f_ptr = c_loc(zarray)
        call fread_zn(comm,f_ptr,nd,ng,nl,nv,path,fileid)
c
        end subroutine fread_z4
c
c       ----------------------------------------------------------------
c
c       Check whether the dataset referred to by path exists
c
c       HDF5 has an subroutine to check existance but if path is more
c       complicated such as group1/group2/dataset then this subroutine
c       will fail if group1 or group2 does not exist. 
c
        logical function fexist(fileid,path)
        implicit none
        integer(hid_t), intent(in) :: fileid ! the file identifier
        character(len=*), intent(in) :: path ! the path to the dataset
        logical :: ofound ! did we find the path?
        integer :: pos ! character position
        integer :: offset ! the last position
        integer(hid_t) :: locid ! locality ID (file or group)
c
        ofound = .true.
        locid = fileid
        offset = 1
        pos = 1
        do while (pos.gt.0 .and. ofound)
          pos = index(path(offset:),"/")
          if (pos.eq.0) then
            call h5lexists_f(locid,path,ofound,h5err)
            if (h5err.ne.0) then
              write(*,*)'ERROR: fexist: h5lexists_f failed:',path
            endif
          else
            call h5lexists_f(locid,path(1:offset+pos-1),
     &                       ofound,h5err)
            if (h5err.ne.0) then
              write(*,*)'ERROR: fexist: h5lexists_f failed:',path
            endif
          endif
          offset = offset+pos
        enddo
        fexist = ofound
        return
        end function fexist
c
c       ----------------------------------------------------------------
c
c       Create the dataset referred to by path
c
c       In HDF5 a path to a data set may have the form
c       group1/group2/.../dataset. If some of the groups in the path
c       do not exist they need to be created before the data set
c       can be created. This subroutine creates all necessary groups
c       and then the dataset (similar to "mkdir -p").
c
        subroutine fcreate(fileid,path,memtype,spaceid,dsetid)
        use hdf5
        implicit none
        integer(hid_t), intent(in) :: fileid ! the file identifier
        character(len=*), intent(in) :: path ! the path to the dataset
        integer(hid_t), intent(in) :: memtype ! data memory type
        integer(hid_t), intent(in) :: spaceid ! file space ID
        integer(hid_t), intent(out) :: dsetid  ! data set ID
c
        logical :: ofound = .true. ! did we find the path?
        integer :: pos ! character position
        integer :: offset ! the last position
        integer(hid_t) :: locid ! locality ID (file or group)
        integer(hid_t) :: grpid ! group ID
        integer(size_t) :: hint ! size hint 
c
        hint = 256
        locid = fileid
        offset = 1
        pos = 1
c
c       The loop terminates when pos.eq.0 (i.e. no more "/" found).
c
        do while (pos.gt.0 .and. ofound)
          pos = index(path(offset:),"/")
          if (pos.eq.0) then
            call h5lexists_f(locid,path,ofound,h5err)
            if (h5err.ne.0) then
              write(*,*)'ERROR: fcreate: h5lexists_f failed:',path
            endif
          else
            call h5lexists_f(locid,path(1:offset+pos-1),
     &                       ofound,h5err)
            if (h5err.ne.0) then
              write(*,*)'ERROR: fcreate: h5lexists_f failed:',path
            endif
          endif
          if (.not. ofound) then
            if (pos.eq.0) then
c
c             no missing groups hence create data set
c
              call h5dcreate_f(fileid,path,memtype,spaceid,
     &                         dsetid,h5err)
              if (h5err.ne.0) then
                write(*,*)'ERROR: fcreate: h5dcreate_f failed:',path
              endif
              ofound = .true.
            else
c
c             missing groups hence create group
c
              call h5gcreate_f(fileid,path(1:offset+pos-1),grpid,
     &                         h5err,hint)
              if (h5err.ne.0) then
                write(*,*)'ERROR: fcreate: h5gcreate_f failed:',path
              endif
              call h5gclose_f(grpid,h5err)
              if (h5err.ne.0) then
                write(*,*)'ERROR: fcreate: h5gclose_f failed:',path
              endif
              ofound = .true.
            endif
          endif
          offset = offset+pos
        enddo
        end subroutine fcreate
c
c       ----------------------------------------------------------------
c
#endif
      end module hdf5io_mod
c
c-----------------------------------------------------------------------
c
#ifdef TEST_HDF5
      program test_hdf5
      use hdf5io_mod
      implicit none
c
c     HDF5 offers parallel I/O capabilities but there is no clear
c     specification of how these operations work. For example
c     if a collective write operation is issued must all processors
c     provide HDF5 with the same data (as in a replicated data
c     parallelization mode) so that HDF5 can pick and choose what
c     part of the data from each processor to save? Or alternatively
c     is the collective operation just saving each processors data
c     and the collective nature of the call simply implies a
c     synchronization? Neither the HDF5 documentation nor the provided
c     code examples provide any information on the semantics of the
c     parallel I/O operations.
c
c     The initial code to support with HDF5 is the GW code by Andrey
c     Kutepov. This code is parallelized using a mix of different
c     models for different variables and code regions. There are
c     3 different parallelization models:
c
c     1. Conventional replicated data parallelism: All processors have
c        their own copy of the full data structure. Write operations
c        are issued by processor 0 storing its copy of the data.
c        Read operations are also issued by processor 0, followed by
c        broadcast to replicated the data to all processors.
c
c     2. Mixed replicated/distributed data parallelism: The pool of
c        processors is broken up into groups. The processors within
c        a group all have identical copies of the relevant data (as in
c        replicated data parallelism). However, large data structures
c        are distributed over groups, so that each group operates on a
c        different part.
c
c     3. Conventional distributed data parallelism: A large data
c        structure is split over all processors, every processor working
c        on its own unique segment of the data. In this case every
c        processor is individually responsible for writing and reading
c        its data segment.
c
c     HDF5 can be successfully introduced only if an approach can be
c     formulated to support all three parallelization strategies. 
c     Ideally the following operations should be possible:
c     A. One data file for all data.
c     B. Write and read replicated data without using additional
c        MPI calls.
c     C. Write and read data replicated in each group without
c        using additional MPI calls.
c     D. Write and read fully distributed data.
c     E. Doing all of the above but writing and reading data on 
c        different number of processors. E.g. preparing the initial
c        data set on 32 processors and then using it on 128 processors.
c
c     The program below performs a number of tests to exercise the
c     required scenarios A-E. The testing approach is to provide
c     write and read routines. Each write invokes HDF5 to store data
c     in a particular way, opening a file, saving the data, and closing
c     the file. Each read invokes HDF5 to retrieve the data, again
c     opening the file, reading the data, and closing the file. 
c     When the read routines return, the data obtained is checked 
c     against what was supposed to be returned. Ensuring that each read
c     and write go through the whole sequence operations from open to
c     close enables invoking different read/write operations in 
c     arbitrary combinations.
c
c     Furthermore, the test data is created in such a way that at the
c     top level routines all processors have complete and identical
c     copies of the data. Various distributed data approaches may be
c     simulated by selecting and copying parts of the whole data set.
c     When a read routine returns the data retrieved can be compared
c     to the appropriate section of the full data set to check that the
c     operation completed as expected.
c
#include "mpif.h"
c
      integer :: ndata  ! number of 8-byte words per processor [Input]
      integer :: wrldsz ! the size of MPI_COMM_WORLD
      integer :: wrldrk ! the rank in MPI_COM_WORLD
      integer :: rowsz  ! the size of a row in a processor grid
      integer :: rowrk  ! the current rank in a row 
      integer :: colsz  ! the size of a column in a processor grid
      integer :: colrk  ! the current rank in a column 
      integer :: nrow   ! the number of processor rows
      integer :: ncol   ! the number of processor columns
      integer :: ii     ! counter
      integer :: jj     ! counter
      integer :: ia     ! character number
      integer :: ip     ! position
      integer :: it     ! index in serialized array
c
      integer :: comm_col ! communicator across rows for column 0
      integer :: comm_row ! communicator across column for each row
      integer :: comm_one ! communicator for each world rank
c
      integer :: mpierr ! error flag
c
      logical :: ok
c
      character(len=32) :: arg
c
      double precision, allocatable :: prim_data(:) ! the primary data
      double precision, allocatable :: read_data(:) ! read back data
      double precision, allocatable :: scr(:) ! scratch
      double complex, allocatable :: prim_cmpl(:) ! the primary data
      double complex, allocatable :: read_cmpl(:) ! read back data
      double complex, allocatable :: scr_cmplx(:) ! scratch
      character*3, allocatable :: prim_char(:,:) ! the primary data
      character*3, allocatable :: read_char(:,:) ! read back data
c
      character(len=52) :: alphabet = 
     &  "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQERSTUVWXYZ"
c
      logical :: prim_l1
      logical :: prim_l2
      logical :: read_l1
      logical :: read_l2
c
      character(len=52) :: c1
      double precision  :: d1
      integer           :: i1
      logical           :: l1
      double complex    :: z1
c
      c1 = alphabet
      d1 = acos(-1.0d0)
      i1 = 222
      l1 = .true.
      z1 = cmplx(d1,asin(-1.0d0))
c
      call mpi_init(mpierr)
      call mpi_comm_size(MPI_COMM_WORLD,wrldsz,mpierr)
      call mpi_comm_rank(MPI_COMM_WORLD,wrldrk,mpierr)
      call find_proc_grid(wrldsz,nrow,ncol)
c
      call finit()
c
c     We assume interconversions between the processor rank and the
c     the position in the processor grid to work the same way as the
c     order in the elements of a matrix.
c
      rowrk = mod(wrldrk,nrow)
      colrk = wrldrk/ncol
      if (wrldrk.eq.0) then
        write(*,*)"world  size = ",wrldsz
        write(*,*)"row    size = ",ncol
        write(*,*)"column size = ",nrow
        if (nrow.eq.wrldsz.or.ncol.eq.wrldsz) then
          write(*,*)"Useless processor grid"
c         call mpi_abort(MPI_COMM_WORLD,mpierr)
        endif
      endif
c
c     Split communicators
c
      call mpi_comm_split(MPI_COMM_WORLD,rowrk,colrk,comm_row,mpierr)
      call mpi_comm_split(MPI_COMM_WORLD,colrk,rowrk,comm_col,mpierr)
      call mpi_comm_split(MPI_COMM_WORLD,wrldrk,wrldrk,comm_one,mpierr)
c     if (colrk.eq.0) then
c       call mpi_comm_split(MPI_COMM_WORLD,colrk,rowrk,comm_col,mpierr)
c     else
c       call mpi_comm_split(MPI_COMM_WORLD,MPI_UNDEFINED,rowrk,comm_col,
c    &                      mpierr)
c       comm_col = MPI_UNDEFINED
c     endif
c
      if (wrldrk.eq.0) then
        call get_command_argument(1,arg)
        read(arg,*)ndata
      endif
      call mpi_bcast(ndata,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
c
      allocate(prim_data(1:ndata))
      allocate(prim_cmpl(1:ndata))
      allocate(prim_char(1:ndata,1:ndata))
      prim_char = "   "
      do ii = 1, ndata
        prim_data(ii) = ii
        prim_cmpl(ii) = complex(ii,ii)
        do jj = 1, ndata
          it = (jj-1)*ndata+ii
          ip = 1+mod(it,3)
          ia = 1+mod(it,52)
          prim_char(ii,jj)(ip:ip) = alphabet(ia:ia)
        enddo
      enddo
      prim_l1 = .true.
      prim_l2 = .false.
      call write_replicated_s("words.dat",MPI_COMM_WORLD,c1,d1,i1,l1,z1)
      call write_replicated("words.dat",MPI_COMM_WORLD,prim_data,ndata)
      call write_replicated2("words.dat",MPI_COMM_WORLD,prim_data,ndata)
      call write_replicated_l("words.dat",MPI_COMM_WORLD,
     &                        prim_l1,prim_l2)
      call write_replicated_t("words.dat",MPI_COMM_WORLD,prim_char,
     &                        ndata)
      call write_distributed_3d("words.dat",comm_one,prim_data,ndata)
      call write_distributed_3z("words.dat",comm_one,prim_cmpl,ndata)
      call rdwrt_replicated_s("words.dat",MPI_COMM_WORLD,c1,d1,i1,l1,z1)
      call read_replicated_s("words.dat",MPI_COMM_WORLD,c1,d1,i1,l1,z1)
      allocate(read_data(1:ndata))
      allocate(read_cmpl(1:ndata))
      allocate(read_char(1:ndata,1:ndata))
      allocate(scr(1:ndata))
      allocate(scr_cmplx(1:ndata))
      read_data=0.0d0
      read_char="   "
      call read_replicated_l("words.dat",MPI_COMM_WORLD,read_l1,read_l2)
      if (wrldrk.eq.0) then
        ok = .true.
        ok = ok.and.(prim_l1.eqv.read_l1)
        ok = ok.and.(prim_l2.eqv.read_l2)
        if (.not.ok) then
          write(*,*)'HvD: repl_l: ',prim_l1,read_l1
          write(*,*)'HvD: repl_l: ',prim_l2,read_l2
        endif
      endif
      call read_replicated_t("words.dat",MPI_COMM_WORLD,read_char,ndata)
      if (wrldrk.eq.0) then
        ok = .true.
        do ii = 1, ndata
          do jj = 1, ndata
            ok = prim_char(ii,jj).eq.read_char(ii,jj)
            if (.not.ok) then
              write(*,*)'HvD: repl_t: ',ii,jj,
     &                  prim_char(ii,jj),read_char(ii,jj)
            endif
          enddo
        enddo
      endif
      call read_replicated("words.dat",MPI_COMM_WORLD,read_data,ndata)
      if (wrldrk.eq.0) then
        ok = .true.
        do ii = 1, ndata
          ok=ok.and.(abs(prim_data(ii)-read_data(ii)).le.1.0d-10)
        enddo
        if (.not.ok) then
          do ii = 1, ndata
            write(*,*)'HvD: repl: ',ii,prim_data(ii),read_data(ii)
          enddo
        endif
      endif
      call mpi_barrier(MPI_COMM_WORLD,mpierr)
      read_data=0.0d0
      scr=0.0d0
      call read_distributed_3d("words.dat",comm_one,
     &                         scr,ndata)
      call mpi_allreduce(scr,read_data,ndata,MPI_DOUBLE,MPI_SUM,
     &                   MPI_COMM_WORLD,mpierr)
      if (wrldrk.eq.0) then
        write(*,*)
        ok = .true.
        do ii = 1, ndata
          ok=ok.and.(abs(prim_data(ii)-read_data(ii)).le.1.0d-10)
        enddo
        if (.not.ok) then
          do ii = 1, ndata
            write(*,*)'HvD: dist: ',ii,prim_data(ii),read_data(ii)
          enddo
        endif
      endif
      call mpi_barrier(MPI_COMM_WORLD,mpierr)
      read_cmpl=complex(0.0d0,0.0d0)
      scr_cmplx=complex(0.0d0,0.0d0)
      call read_distributed_3z("words.dat",comm_one,
     &                         scr_cmplx,ndata)
      call mpi_allreduce(scr_cmplx,read_cmpl,ndata,MPI_DOUBLE_COMPLEX,
     &                   MPI_SUM,MPI_COMM_WORLD,mpierr)
      if (wrldrk.eq.0) then
        write(*,*)
        ok = .true.
        do ii = 1, ndata
          ok=ok.and.(abs(prim_cmpl(ii)-read_cmpl(ii)).le.1.0d-10)
        enddo
        if (.not.ok) then
          do ii = 1, ndata
            write(*,*)'HvD: dist: ',ii,prim_cmpl(ii),read_cmpl(ii)
          enddo
        endif
      endif
      call mpi_barrier(MPI_COMM_WORLD,mpierr)
      read_data=0.0d0
c     call read_replibuted("words.dat",comm_row,read_data,ndata)
c     if (wrldrk.eq.0) then
c       write(*,*)
c       write(*,*)'HvD: prim_data=',prim_data
c     endif
c     flush(6)
      call mpi_barrier(MPI_COMM_WORLD,mpierr)
c     flush(6)
c     write(*,*)wrldrk,' HvD: read_data=',read_data
c TESTS HERE
      deallocate(read_char)
      deallocate(prim_char)
      deallocate(scr_cmplx)
      deallocate(read_cmpl)
      deallocate(prim_cmpl)
      deallocate(scr)
      deallocate(read_data)
      deallocate(prim_data)
c
      call ffinal()
c
      call mpi_comm_free(comm_one,mpierr)
      call mpi_comm_free(comm_row,mpierr)
      call mpi_comm_free(comm_col,mpierr)
c
      call mpi_finalize(mpierr)
c
      CONTAINS
c
c       ----------------------------------------------------------------
c
        subroutine write_replicated_s(filename,comm,c1,d1,i1,l1,z1)
        use hdf5
        use hdf5io_mod
        implicit none
#include "mpif.h"
c
c       The data passed are replicated scalar data structures. I.e.
c       every rank holds the same data. The file is opened or created
c       the data is written and the file is closed.
c
        character*(*), intent(in) :: filename
        integer, intent(in) :: comm        ! the communicator to use
        character*(*),    intent(in) :: c1 ! a string scalar
        double precision, intent(in) :: d1 ! a double scalar
        integer,          intent(in) :: i1 ! an integer scalar
        logical,          intent(in) :: l1 ! a boolean scalar
        double complex,   intent(in) :: z1 ! a complex scalar
c
c       Local
c
        integer :: wrld_size
        integer :: wrld_rank
        integer :: mpierr
c
        integer(hid_t) :: fileid
c
c       Code
c
        call fopen(filename,fileid)
        call fwrite(comm,c1,"char1",fileid)
        call fwrite(comm,d1,"double1",fileid)
        call fwrite(comm,i1,"integer1",fileid)
        call fwrite(comm,l1,"logical1",fileid)
        call fwrite(comm,z1,"complex1",fileid)
        call fclose(fileid)
c
        end subroutine write_replicated_s
c
c       ----------------------------------------------------------------
c
        subroutine read_replicated_s(filename,comm,c1,d1,i1,l1,z1)
        use hdf5
        use hdf5io_mod
        implicit none
#include "mpif.h"
c
c       The data passed are replicated scalar data structures. I.e.
c       every rank holds the same data. The file is opened or created
c       the data is written and the file is closed.
c
        character*(*), intent(in) :: filename
        integer, intent(in) :: comm        ! the communicator to use
        character*(*),    intent(in) :: c1 ! a string scalar
        double precision, intent(in) :: d1 ! a double scalar
        integer,          intent(in) :: i1 ! an integer scalar
        logical,          intent(in) :: l1 ! a boolean scalar
        double complex,   intent(in) :: z1 ! a complex scalar
c
c       Local
c
        integer :: wrld_size
        integer :: wrld_rank
        integer :: mpierr
c
        character*(52)   :: cl ! a string scalar
        double precision :: dl ! a double scalar
        integer          :: il ! an integer scalar
        logical          :: ll ! a boolean scalar
        double complex   :: zl ! a complex scalar
c
        integer(hid_t) :: fileid
c
c       Code
c
        call fopen(filename,fileid)
        call fread(comm,cl,"char1",fileid)
        call fread(comm,dl,"double1",fileid)
        call fread(comm,il,"integer1",fileid)
        call fread(comm,ll,"logical1",fileid)
        call fread(comm,zl,"complex1",fileid)
        call fclose(fileid)
c
        if (cl.ne.c1) then
          write(*,*)'cl != c1',cl,c1
        endif
        if (dl.ne.d1) then
          write(*,*)'dl != d1',dl,d1
        endif
        if (il.ne.i1) then
          write(*,*)'il != i1',il,i1
        endif
        if (ll.neqv.l1) then
          write(*,*)'ll != l1',ll,l1
        endif
        if (zl.ne.z1) then
          write(*,*)'zl != z1',zl,z1
        endif
c
        end subroutine read_replicated_s
c
c       ----------------------------------------------------------------
c
        subroutine rdwrt_replicated_s(filename,comm,c1,d1,i1,l1,z1)
        use hdf5
        use hdf5io_mod
        implicit none
#include "mpif.h"
c
c       The data passed are replicated scalar data structures. I.e.
c       every rank holds the same data. The file is opened or created
c       the data is written and the file is closed.
c
        character*(*), intent(in) :: filename
        integer, intent(in) :: comm        ! the communicator to use
        character*(*),    intent(in) :: c1 ! a string scalar
        double precision, intent(in) :: d1 ! a double scalar
        integer,          intent(in) :: i1 ! an integer scalar
        logical,          intent(in) :: l1 ! a boolean scalar
        double complex,   intent(in) :: z1 ! a complex scalar
c
c       Local
c
        integer :: wrld_size
        integer :: wrld_rank
        integer :: mpierr
c
        character*(52)   :: cl ! a string scalar
        double precision :: dl ! a double scalar
        integer          :: il ! an integer scalar
        logical          :: ll ! a boolean scalar
        double complex   :: zl ! a complex scalar
c
        integer(hid_t) :: fileid
c
c       Code
c
        call fopen(filename,fileid)
        call fwrite(comm,c1,"chara",fileid)
        call fwrite(comm,d1,"doublea",fileid)
        call fwrite(comm,i1,"integera",fileid)
        call fwrite(comm,l1,"logicala",fileid)
        call fwrite(comm,z1,"complexa",fileid)
        call fread(comm,cl,"chara",fileid)
        call fread(comm,dl,"doublea",fileid)
        call fread(comm,il,"integera",fileid)
        call fread(comm,ll,"logicala",fileid)
        call fread(comm,zl,"complexa",fileid)
        call fclose(fileid)
c
        if (cl.ne.c1) then
          write(*,*)'cl != c1',cl,c1
        endif
        if (dl.ne.d1) then
          write(*,*)'dl != d1',dl,d1
        endif
        if (il.ne.i1) then
          write(*,*)'il != i1',il,i1
        endif
        if (ll.neqv.l1) then
          write(*,*)'ll != l1',ll,l1
        endif
        if (zl.ne.z1) then
          write(*,*)'zl != z1',zl,z1
        endif
c
        end subroutine rdwrt_replicated_s
c
c       ----------------------------------------------------------------
c
        subroutine write_distributed_3d(filename,comm,words,nw)
        use hdf5
        use hdf5io_mod
        implicit none
#include "mpif.h"
c
c       The data passed is a replicated data structure. That means every
c       rank presents exactly the same data. Also the filename is given.
c       Hence we create the HDF5 file, collectively write the data, and
c       close the file.
c
        character*(*), intent(in) :: filename
        integer, intent(in) :: comm ! the communicator to use
        integer, intent(in) :: nw   ! the number of double precision words
        double precision, intent(in), target :: words(nw) ! the words of data
c
c       Local
c
        integer :: wrld_size
        integer :: wrld_rank
        integer :: nsize
        integer :: irank
        integer :: mpierr
        integer :: dimw ! words per proc
        integer :: dimr ! no. remainder procs
        integer :: dimt ! temporary: 1 if this proc get 1 extra, 0 else
        integer :: offs
        integer :: ng(3)
        integer :: nv(3)
        double precision, allocatable :: rdata(:,:,:)
        double precision, pointer :: rwords(:,:,:)
c
        integer(hid_t) :: fileid
c
c       Code
c
        call mpi_comm_size(MPI_COMM_WORLD,wrld_size,mpierr)
        call mpi_comm_rank(MPI_COMM_WORLD,wrld_rank,mpierr)
        call mpi_comm_size(comm,nsize,mpierr)
        call mpi_comm_rank(comm,irank,mpierr)
c
        dimw = nw/(2*2)/wrld_size
        dimr = mod(nw/(2*2),wrld_size)
        dimt = 0
        if (wrld_rank.lt.dimr) dimt=1
        if (dimw+dimt.le.0) then
          write(*,*)'write_distributed_3d: num words too small'
          write(*,*)'num words     = ',nw
          write(*,*)'num MPI ranks = ',wrld_size
          write(*,*)'MPI rank      = ',wrld_rank
          write(*,*)'two more dims = ',2,2
c         call mpi_abort(MPI_COMM_WORLD,mpierr)
        endif
        if (mod(nw,2*2).ne.0) then
          write(*,*)'write_distributed_3d: not a multiple'
          write(*,*)'num words     = ',nw
          write(*,*)'num MPI ranks = ',wrld_size
          write(*,*)'two more dims = ',2,2
          call mpi_abort(MPI_COMM_WORLD,mpierr)
        endif
        allocate(rdata(2,dimw+dimt,2))
        rwords(1:2,1:nw/(2*2),1:2)=>words
        if (wrld_rank.lt.dimr) then
          rdata(:,1:dimw+dimt,:)=rwords(:,(dimw+dimt)*wrld_rank+1:
     &                                    (dimw+dimt)*(wrld_rank+1),:)
        else
          rdata(:,1:dimw+dimt,:)=rwords(:,
     &         (dimw+1)*dimr+dimw*(wrld_rank-dimr)+1:
     &         (dimw+1)*dimr+dimw*(wrld_rank-dimr+1),:)
        endif
c
        call fopen(filename,fileid)
        nv = 0
        nv(2) = dimw*wrld_rank
        if (wrld_rank.lt.dimr) nv(2)=nv(2)+wrld_rank
        ng(1) = 2
        ng(2) = nw/(2*2)
        ng(3) = 2
        call fwrite(comm,rdata,ng,nv,"rdata1",fileid)
        call fclose(fileid)
c
        deallocate(rdata)
c
        end subroutine write_distributed_3d
c
c       ----------------------------------------------------------------
c
        subroutine read_distributed_3d(filename,comm,words,nw)
        use hdf5
        use hdf5io_mod
        implicit none
#include "mpif.h"
c
c       The data passed is a replicated data structure. That means every
c       rank presents exactly the same data. Also the filename is given.
c       Hence we create the HDF5 file, collectively write the data, and
c       close the file.
c
        character*(*), intent(in) :: filename
        integer, intent(in) :: comm ! the communicator to use
        integer, intent(in) :: nw   ! the number of double precision words
        double precision, intent(out), target :: words(nw) ! the words of data
c
c       Local
c
        integer :: wrld_size
        integer :: wrld_rank
        integer :: nsize
        integer :: irank
        integer :: mpierr
        integer :: dimw ! words per proc
        integer :: dimr ! no. remainder procs
        integer :: dimt ! temporary: 1 if this proc get 1 extra, 0 else
        integer :: offs
        integer :: ng(3)
        integer :: nv(3)
        double precision, allocatable :: rdata(:,:,:)
        double precision, pointer :: rwords(:,:,:)
c
        integer(hid_t) :: fileid
c
c       Code
c
        call mpi_comm_size(MPI_COMM_WORLD,wrld_size,mpierr)
        call mpi_comm_rank(MPI_COMM_WORLD,wrld_rank,mpierr)
        call mpi_comm_size(comm,nsize,mpierr)
        call mpi_comm_rank(comm,irank,mpierr)
c
        dimw = nw/(2*2)/wrld_size
        dimr = mod(nw/(2*2),wrld_size)
        dimt = 0
        if (wrld_rank.lt.dimr) dimt=1
        if (dimw+dimt.le.0) then
          write(*,*)'read_distributed_3d: num words too small'
          write(*,*)'num words     = ',nw
          write(*,*)'num MPI ranks = ',wrld_size
          write(*,*)'MPI ranks     = ',wrld_rank
          write(*,*)'two more dims = ',2,2
c         call mpi_abort(MPI_COMM_WORLD,mpierr)
        endif
        if (mod(nw,2*2).ne.0) then
          write(*,*)'read_distributed_3d: not a multiple'
          write(*,*)'num words     = ',nw
          write(*,*)'num MPI ranks = ',wrld_size
          write(*,*)'two more dims = ',2,2
          call mpi_abort(MPI_COMM_WORLD,mpierr)
        endif
        allocate(rdata(2,dimw+dimt,2))
        rwords(1:2,1:nw/(2*2),1:2)=>words
c
        call fopen(filename,fileid)
        nv = 0
        nv(2) = dimw*wrld_rank
        if (wrld_rank.lt.dimr) nv(2)=nv(2)+wrld_rank
        ng(1) = 2
        ng(2) = nw/(2*2)
        ng(3) = 2
        call fread(comm,rdata,ng,nv,"rdata1",fileid)
        call fclose(fileid)
        if (wrld_rank.lt.dimr) then
          rwords(:,(dimw+dimt)*wrld_rank+1:(dimw+dimt)*(wrld_rank+1),:)
     &    =rdata(:,1:dimw+dimt,:)
        else
          rwords(:,(dimw+1)*dimr+dimw*(wrld_rank-dimr)+1:
     &             (dimw+1)*dimr+dimw*(wrld_rank-dimr+1),:)
     &    =rdata(:,1:dimw+dimt,:)
        endif
c
        deallocate(rdata)
c
        end subroutine read_distributed_3d
c
c       ----------------------------------------------------------------
c
        subroutine write_distributed_3z(filename,comm,words,nw)
        use hdf5
        use hdf5io_mod
        implicit none
#include "mpif.h"
c
c       The data passed is a replicated data structure. That means every
c       rank presents exactly the same data. Also the filename is given.
c       Hence we create the HDF5 file, collectively write the data, and
c       close the file.
c
        character*(*), intent(in) :: filename
        integer, intent(in) :: comm ! the communicator to use
        integer, intent(in) :: nw   ! the number of double precision words
        double complex, intent(in), target :: words(nw) ! the words of data
c
c       Local
c
        integer :: wrld_size
        integer :: wrld_rank
        integer :: nsize
        integer :: irank
        integer :: mpierr
        integer :: dimw ! words per proc
        integer :: dimr ! no. remainder procs
        integer :: dimt ! temporary: 1 if this proc get 1 extra, 0 else
        integer :: offs
        integer :: ng(3)
        integer :: nv(3)
        double complex, allocatable :: rdata(:,:,:)
        double complex, pointer :: rwords(:,:,:)
c
        integer(hid_t) :: fileid
c
c       Code
c
        call mpi_comm_size(MPI_COMM_WORLD,wrld_size,mpierr)
        call mpi_comm_rank(MPI_COMM_WORLD,wrld_rank,mpierr)
        call mpi_comm_size(comm,nsize,mpierr)
        call mpi_comm_rank(comm,irank,mpierr)
c
        dimw = nw/(2*2)/wrld_size
        dimr = mod(nw/(2*2),wrld_size)
        dimt = 0
        if (wrld_rank.lt.dimr) dimt=1
        if (dimw+dimt.le.0) then
          write(*,*)'write_distributed_3z: num words too small'
          write(*,*)'num words     = ',nw
          write(*,*)'num MPI ranks = ',wrld_size
          write(*,*)'MPI rank      = ',wrld_rank
          write(*,*)'two more dims = ',2,2
c         call mpi_abort(MPI_COMM_WORLD,mpierr)
        endif
        if (mod(nw,2*2).ne.0) then
          write(*,*)'write_distributed_3z: not a multiple'
          write(*,*)'num words     = ',nw
          write(*,*)'num MPI ranks = ',wrld_size
          write(*,*)'two more dims = ',2,2
          call mpi_abort(MPI_COMM_WORLD,mpierr)
        endif
        allocate(rdata(2,dimw+dimt,2))
        rwords(1:2,1:nw/(2*2),1:2)=>words
        if (wrld_rank.lt.dimr) then
          rdata(:,1:dimw+dimt,:)=rwords(:,(dimw+dimt)*wrld_rank+1:
     &                                    (dimw+dimt)*(wrld_rank+1),:)
        else
          rdata(:,1:dimw+dimt,:)=rwords(:,
     &         (dimw+1)*dimr+dimw*(wrld_rank-dimr)+1:
     &         (dimw+1)*dimr+dimw*(wrld_rank-dimr+1),:)
        endif
c
        call fopen(filename,fileid)
        nv = 0
        nv(2) = dimw*wrld_rank
        if (wrld_rank.lt.dimr) nv(2)=nv(2)+wrld_rank
        ng(1) = 2
        ng(2) = nw/(2*2)
        ng(3) = 2
        call fwrite(comm,rdata,ng,nv,"zdata1",fileid)
        call fclose(fileid)
c
        deallocate(rdata)
c
        end subroutine write_distributed_3z
c
c       ----------------------------------------------------------------
c
        subroutine read_distributed_3z(filename,comm,words,nw)
        use hdf5
        use hdf5io_mod
        implicit none
#include "mpif.h"
c
c       The data passed is a replicated data structure. That means every
c       rank presents exactly the same data. Also the filename is given.
c       Hence we create the HDF5 file, collectively write the data, and
c       close the file.
c
        character*(*), intent(in) :: filename
        integer, intent(in) :: comm ! the communicator to use
        integer, intent(in) :: nw   ! the number of double precision words
        double complex, intent(out), target :: words(nw) ! the words of data
c
c       Local
c
        integer :: wrld_size
        integer :: wrld_rank
        integer :: nsize
        integer :: irank
        integer :: mpierr
        integer :: dimw ! words per proc
        integer :: dimr ! no. remainder procs
        integer :: dimt ! temporary: 1 if this proc get 1 extra, 0 else
        integer :: offs
        integer :: ng(3)
        integer :: nv(3)
        double complex, allocatable :: rdata(:,:,:)
        double complex, pointer :: rwords(:,:,:)
c
        integer(hid_t) :: fileid
c
c       Code
c
        call mpi_comm_size(MPI_COMM_WORLD,wrld_size,mpierr)
        call mpi_comm_rank(MPI_COMM_WORLD,wrld_rank,mpierr)
        call mpi_comm_size(comm,nsize,mpierr)
        call mpi_comm_rank(comm,irank,mpierr)
c
        dimw = nw/(2*2)/wrld_size
        dimr = mod(nw/(2*2),wrld_size)
        dimt = 0
        if (wrld_rank.lt.dimr) dimt=1
        if (dimw+dimt.le.0) then
          write(*,*)'read_distributed_3z: num words too small'
          write(*,*)'num words     = ',nw
          write(*,*)'num MPI ranks = ',wrld_size
          write(*,*)'MPI ranks     = ',wrld_rank
          write(*,*)'two more dims = ',2,2
c         call mpi_abort(MPI_COMM_WORLD,mpierr)
        endif
        if (mod(nw,2*2).ne.0) then
          write(*,*)'read_distributed_3z: not a multiple'
          write(*,*)'num words     = ',nw
          write(*,*)'num MPI ranks = ',wrld_size
          write(*,*)'two more dims = ',2,2
          call mpi_abort(MPI_COMM_WORLD,mpierr)
        endif
        allocate(rdata(2,dimw+dimt,2))
        rwords(1:2,1:nw/(2*2),1:2)=>words
c
        call fopen(filename,fileid)
        nv = 0
        nv(2) = dimw*wrld_rank
        if (wrld_rank.lt.dimr) nv(2)=nv(2)+wrld_rank
        ng(1) = 2
        ng(2) = nw/(2*2)
        ng(3) = 2
        call fread(comm,rdata,ng,nv,"zdata1",fileid)
        call fclose(fileid)
        if (wrld_rank.lt.dimr) then
          rwords(:,(dimw+dimt)*wrld_rank+1:(dimw+dimt)*(wrld_rank+1),:)
     &    =rdata(:,1:dimw+dimt,:)
        else
          rwords(:,(dimw+1)*dimr+dimw*(wrld_rank-dimr)+1:
     &             (dimw+1)*dimr+dimw*(wrld_rank-dimr+1),:)
     &    =rdata(:,1:dimw+dimt,:)
        endif
c
        deallocate(rdata)
c
        end subroutine read_distributed_3z
c
c       ----------------------------------------------------------------
c
      end program test_hdf5
c
c-----------------------------------------------------------------------
c
      subroutine write_replicated(filename,comm,words,nw)
      use hdf5
      use hdf5io_mod
      implicit none
#include "mpif.h"
c
c     The data passed is a replicated data structure. That means every
c     rank presents exactly the same data. Also the filename is given.
c     Hence we create the HDF5 file, collectively write the data, and
c     close the file.
c
      character*(*), intent(in) :: filename
      integer, intent(in) :: comm ! the communicator to use
      integer, intent(in) :: nw   ! the number of double precision words
      double precision, intent(in) :: words(nw) ! the words of data
c
c     Local
c
      integer :: wrld_nsize
      integer :: wrld_irank
      integer :: nsize
      integer :: irank
      integer :: mpierr
      integer :: dimw
      integer :: offs
c
      integer(hid_t) :: fileid
c
c     Code
c
      call mpi_comm_size(MPI_COMM_WORLD,wrld_nsize,mpierr)
      call mpi_comm_rank(MPI_COMM_WORLD,wrld_irank,mpierr)
      call mpi_comm_size(comm,nsize,mpierr)
      call mpi_comm_rank(comm,irank,mpierr)
c
      dimw = nw
      offs = 0
c
      call fopen(filename,fileid)
      call fwrite(comm,words(offs+1:offs+dimw),nw,offs,"data1",fileid)
      call fclose(fileid)
c
      end subroutine write_replicated
c
c-----------------------------------------------------------------------
c
      subroutine write_replicated_l(filename,comm,log_l1,log_l2)
      use hdf5
      use hdf5io_mod
      implicit none
#include "mpif.h"
c
c     The data passed is a replicated data structure. That means every
c     rank presents exactly the same data. Also the filename is given.
c     Hence we create the HDF5 file, collectively write the data, and
c     close the file.
c
      character*(*), intent(in) :: filename
      integer, intent(in) :: comm ! the communicator to use
      logical, intent(in) :: log_l1 ! logical value number 1
      logical, intent(in) :: log_l2 ! logical value number 2
c
c     Local
c
      integer :: wrld_nsize
      integer :: wrld_irank
      integer :: nsize
      integer :: irank
      integer :: mpierr
c
      integer(hid_t) :: fileid
c
c     Code
c
      call mpi_comm_size(MPI_COMM_WORLD,wrld_nsize,mpierr)
      call mpi_comm_rank(MPI_COMM_WORLD,wrld_irank,mpierr)
      call mpi_comm_size(comm,nsize,mpierr)
      call mpi_comm_rank(comm,irank,mpierr)
c
      call fopen(filename,fileid)
      call fwrite(comm,log_l1,"log_l1",fileid)
      call fwrite(comm,log_l2,"log_l2",fileid)
      call fclose(fileid)
c
      end subroutine write_replicated_l
c
c-----------------------------------------------------------------------
c
      subroutine write_replicated_t(filename,comm,text,ndata)
      use hdf5
      use hdf5io_mod
      implicit none
#include "mpif.h"
c
c     The data passed is a replicated data structure. That means every
c     rank presents exactly the same data. Also the filename is given.
c     Hence we create the HDF5 file, collectively write the data, and
c     close the file.
c
      character*(*), intent(in) :: filename
      integer, intent(in) :: comm ! the communicator to use
      integer, intent(in) :: ndata ! the dimension
      character(len=*), intent(in) :: text(ndata,ndata) ! characters
c
c     Local
c
      integer :: wrld_nsize
      integer :: wrld_irank
      integer :: nsize
      integer :: irank
      integer :: mpierr
c
      integer(hid_t) :: fileid
      integer :: nd
      integer :: ng(2)
      integer :: nv(2)
c
c     Code
c
      call mpi_comm_size(MPI_COMM_WORLD,wrld_nsize,mpierr)
      call mpi_comm_rank(MPI_COMM_WORLD,wrld_irank,mpierr)
      call mpi_comm_size(comm,nsize,mpierr)
      call mpi_comm_rank(comm,irank,mpierr)
c
      ng(1) = size(text,1)
      ng(2) = size(text,2)
      nv    = 0
c
      call fopen(filename,fileid)
      call fwrite(comm,text,ng,nv,"text",fileid)
      call fclose(fileid)
c
      end subroutine write_replicated_t
c
c-----------------------------------------------------------------------
c
      subroutine write_replicated2(filename,comm,words,nw)
      use hdf5
      use hdf5io_mod
      implicit none
c
c     The data passed is a replicated data structure. That means every
c     rank presents exactly the same data. Also the filename is given.
c     Hence we create the HDF5 file, collectively write the data, and
c     close the file.
c
      character*(*), intent(in) :: filename
      integer, intent(in) :: comm ! the communicator to use
      integer, intent(in) :: nw   ! the number of double precision words
      double precision, intent(in) :: words(nw) ! the words of data
c
c     Local
c
      integer :: nsize
      integer :: irank
      integer :: mpierr
      integer :: dimw
      integer :: offs
      integer :: idata
c
      integer(hid_t) :: fileid
c
c     Code
c
      call mpi_comm_size(comm,nsize,mpierr)
      call mpi_comm_rank(comm,irank,mpierr)
c
      idata = words(1)
c
      call fopen(filename,fileid)
      call fwrite(comm,idata,"idata1",fileid)
      call fclose(fileid)
c
      end subroutine write_replicated2
c
c-----------------------------------------------------------------------
c
      subroutine read_replicated(filename,comm,words,nw)
      use hdf5
      use hdf5io_mod
      implicit none
c
c     The data passed is a replicated data structure. That means every
c     rank returns exactly the same data. Also the filename is given.
c     Hence we open the HDF5 file, collectively read the data, and
c     close the file.
c
      character*(*), intent(in) :: filename
      integer, intent(in) :: comm ! the communicator to use
      integer, intent(in) :: nw   ! the number of double precision words
      double precision, intent(out) :: words(nw) ! the words of data
c
c     Local
c
      integer :: nsize
      integer :: irank
      integer :: mpierr
      integer :: dimw
      integer :: offs
c
      integer(hid_t) :: fileid
c
c     Code
c
      call mpi_comm_size(comm,nsize,mpierr)
      call mpi_comm_rank(comm,irank,mpierr)
c
      dimw = nw
      offs = 0
c
      call fopen(filename,fileid)
      call fread(comm,words(offs+1:offs+dimw),nw,offs,"data1",fileid)
      call fclose(fileid)
c
      end subroutine read_replicated
c
c-----------------------------------------------------------------------
c
      subroutine read_replicated_l(filename,comm,log_l1,log_l2)
      use hdf5
      use hdf5io_mod
      implicit none
c
c     The data passed is a replicated data structure. That means every
c     rank returns exactly the same data. Also the filename is given.
c     Hence we open the HDF5 file, collectively read the data, and
c     close the file.
c
      character*(*), intent(in) :: filename
      integer, intent(in) :: comm ! the communicator to use
      logical, intent(out) :: log_l1 ! logical value number 1
      logical, intent(out) :: log_l2 ! logical value number 2
c
c     Local
c
      integer :: nsize
      integer :: irank
      integer :: mpierr
c
      integer(hid_t) :: fileid
c
c     Code
c
      call mpi_comm_size(comm,nsize,mpierr)
      call mpi_comm_rank(comm,irank,mpierr)
c
      call fopen(filename,fileid)
      call fread(comm,log_l2,"log_l2",fileid)
      call fread(comm,log_l1,"log_l1",fileid)
      call fclose(fileid)
c
      end subroutine read_replicated_l
c
c-----------------------------------------------------------------------
c
      subroutine read_replicated_t(filename,comm,read_char,ndata)
      use hdf5
      use hdf5io_mod
      implicit none
c
c     The data passed is a replicated data structure. That means every
c     rank returns exactly the same data. Also the filename is given.
c     Hence we open the HDF5 file, collectively read the data, and
c     close the file.
c
      character*(*), intent(in) :: filename
      integer, intent(in) :: comm ! the communicator to use
      integer, intent(in) :: ndata ! the dimensions
      character(len=*), intent(out) :: read_char(ndata,ndata) ! text
c
c     Local
c
      integer :: nsize
      integer :: irank
      integer :: mpierr
      integer :: ng(2)
      integer :: nv(2)
c
      integer(hid_t) :: fileid
c
c     Code
c
      call mpi_comm_size(comm,nsize,mpierr)
      call mpi_comm_rank(comm,irank,mpierr)
c
      ng = ndata
      nv = 0
      call fopen(filename,fileid)
      call fread(comm,read_char,ng,nv,"text",fileid)
      call fclose(fileid)
c
      end subroutine read_replicated_t
c
c-----------------------------------------------------------------------
c
      subroutine read_distributed(filename,comm,words,nw)
      use hdf5
      use hdf5io_mod
      implicit none
c
c     The data passed is a replicated data structure. However, the read
c     operates in a distributed data way. Hence every rank will return
c     only its part of the data.
c
      character*(*), intent(in) :: filename
      integer, intent(in) :: comm ! the communicator to use
      integer, intent(in) :: nw   ! the number of double precision words
      double precision, intent(out) :: words(nw) ! the words of data
c
c     Local
c
      integer :: nsize
      integer :: irank
      integer :: mpierr
      integer :: dimw
      integer :: offs
c
      integer(hid_t) :: fileid
c
c     Code
c
      call mpi_comm_size(comm,nsize,mpierr)
      call mpi_comm_rank(comm,irank,mpierr)
c
      dimw = (nw+nsize-1)/nsize
      offs = irank*dimw
c
      call fopen(filename,fileid)
      call fread(comm,words(offs+1:offs+dimw),nw,offs,"data1",fileid)
      call fclose(fileid)
c
      end subroutine read_distributed
c
c-----------------------------------------------------------------------
c
      subroutine read_replibuted(filename,comm,words,nw)
      use hdf5
      use hdf5io_mod
      implicit none
      include "mpif.h"
c
c     The data passed is a replicated data structure. However in this
c     context the read operation treats the data as partially 
c     replicated and partially distributed. Hence each rank will return
c     the part it replicates, but not the part that is distributed to
c     other ranks.
c
      character*(*), intent(in) :: filename
      integer, intent(in) :: comm ! the communicator to use
      integer, intent(in) :: nw   ! the number of double precision words
      double precision, intent(out) :: words(nw) ! the words of data
c
c     Local
c
      integer :: nsize  ! the size of the comm communicator
      integer :: irank  ! the rank in the comm communicator
      integer :: wrldsz ! the size of the MPI_COMM_WORLD communicator
      integer :: wrldrk ! the rank in the MPI_COMM_WORLD communicator
      integer :: psize  ! the number of processor groups
      integer :: prank  ! the rank of the processor group my rank 
                        ! belongs to
      integer :: mpierr
      integer :: dimw
      integer :: offs
c
      integer(hid_t) :: fileid
c
c     Code
c
      call mpi_comm_size(MPI_COMM_WORLD,wrldsz,mpierr)
      call mpi_comm_rank(MPI_COMM_WORLD,wrldrk,mpierr)
      call mpi_comm_size(comm,nsize,mpierr)
      call mpi_comm_rank(comm,irank,mpierr)
c
      psize = wrldsz/nsize
      prank = wrldrk/psize
c
      dimw = (nw+psize-1)/psize
      offs = prank*dimw
c
      call fopen(filename,fileid)
      call fread(comm,words(offs+1:offs+dimw),nw,offs,"data1",fileid)
      call fclose(fileid)
c
      end subroutine read_replibuted
c
c-----------------------------------------------------------------------
c
      subroutine find_proc_grid(nproc,nrow,ncol)
      implicit none
      integer, intent(in)  :: nproc ! the total number of processors
      integer, intent(out) :: nrow  ! the number of rows
      integer, intent(out) :: ncol  ! the number of columns
c
c     Factorize the number of processors into a number of rows and a
c     number of columns to define a processor grid. Try to obtain a 
c     grid that is as close to square as possible.
c
      integer :: itop ! the top of the loop
      integer :: ii   ! the loop counter
c
      itop = sqrt(1.0d0*nproc)
      itop = min(itop,nproc)
      do ii = itop, 1, -1
        nrow = nproc/ii
        ncol = nproc/nrow
        if (nrow*ncol.eq.nproc) exit
      enddo
c
      end subroutine find_proc_grid
c
c-----------------------------------------------------------------------
#endif
