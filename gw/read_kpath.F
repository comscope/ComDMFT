      subroutine read_kpath(mx_symb,mx_k,n_symb,n_k,s,v_symb,i_symb)
      use manager_mod
      use parallel_mod
      use inp_mod
      implicit none
#ifdef MPI
      include "mpif.h"
#endif
c
c     If the file "kpath" exists use its contents to initialize the
c     band structure k-path. If the file does not exist return without
c     changing any of the arguments.
c
c     In parallel calculations only the root process reads the k-path
c     file. 
c
c     An example of a k-path file is (numbers truncated for space
c     limits):
c
c       X   0.00  -0.50   0.00       G   0.00   0.00   0.00
c       Y   0.50   0.00   0.00       L   0.50  -0.50   0.00
c       G   0.00   0.00   0.00       Z  -0.50   0.00   0.50
c       N  -0.50  -0.50   0.50       G   0.00   0.00   0.00
c       M   0.00   0.00   0.50       R   0.00  -0.50   0.50
c       G   0.00   0.00   0.00
c
      integer, intent(in) ::  mx_symb ! the max number of high symmetry
                                      ! points
      integer, intent(in) ::  mx_k    ! the max number of k-points on
                                      ! the path
      integer, intent(inout) :: n_symb ! the number of high symmetry 
                                       ! points
      integer, intent(inout) :: n_k    ! the number of k-points on the
                                       ! path
      character*3, intent(inout) :: s(mx_symb) ! the high symmetry
                                               ! point labels
      double precision, intent(inout) :: v_symb(3,mx_symb) ! the k-point
                                                           ! coordinates
      integer, intent(inout) :: i_symb(mx_k) ! the k-point indexes on
                                             ! the path
c
      logical :: kpath_exists
      integer :: nk, ik, k_symb ! counters
      integer :: ierr           ! error code
      character*3      :: s1    ! high symmetry point label as read
      double precision :: v1(3) ! k-point coordinates as read
      double precision :: tol = 2.0d-8 ! tolerance 
c
      nk=0
      kpath_exists = .false.
      if (maswrk) then
        inquire(file="kpath",exist=kpath_exists)
      endif
#ifdef MPI
      if (goparr) then
        call brdcst(msgint,kpath_exists,4,master,MPI_COMM_WORLD)
      endif
#endif
      if (.not.kpath_exists) return
c
      if (maswrk) then
        open(3,file="kpath",status="old",form="formatted")
        n_symb = 0
        n_k    = 0
        do while (.true.)
          call inp_line(3,ierr)
          if (ierr.eq.ERR_EOF) exit
          if      (inp_len().eq.0) then
            nk = 0 ! allow for empty lines
          else if (inp_len().eq.4) then
            nk = 1
          else if (inp_len().eq.8) then
            nk = 2
          else
            write(iun,*)'ERROR: read_kpath: inp_line messed up'
          endif
          do ik = 1, nk
            call inp_str(s1,ierr)
            call inp_dbl(v1(1),ierr)
            call inp_dbl(v1(2),ierr)
            call inp_dbl(v1(3),ierr)
            do k_symb = 1, n_symb
              if (s(k_symb).eq.s1.and.
     &            abs(v_symb(1,k_symb)-v1(1)).lt.tol.and.
     &            abs(v_symb(2,k_symb)-v1(2)).lt.tol.and.
     &            abs(v_symb(3,k_symb)-v1(3)).lt.tol) then
c
c               We have this high symmetry point already
c
                n_k = n_k + 1
                if (n_k.gt.mx_k) then
                  write(iun,*)'ERROR: read_kpath: mx_k exceeded'
                  flush(iun)
                  call ending
                endif
                i_symb(n_k) = k_symb
                exit
              endif
            enddo
            if (k_symb.gt.n_symb) then
c
c             We do not have this high symmetry point yet
c
              n_symb = n_symb + 1
              n_k    = n_k + 1
              if (n_symb.gt.mx_symb) then
                write(iun,*)'ERROR: read_kpath: mx_symb exceeded'
                flush(iun)
                call ending
              endif
              if (n_k.gt.mx_k) then
                write(iun,*)'ERROR: read_kpath: mx_k exceeded'
                flush(iun)
                call ending
              endif
              s(n_symb)        = s1
              v_symb(:,n_symb) = v1
              i_symb(n_k)      = n_symb
            endif
          enddo
        enddo
        close(3)
        n_k = n_k - 1 ! The no. paths is one less than the no.
                      ! of high symmetry points (every path 
                      ! connects two consecutive points)
      endif
#ifdef MPI
      if (goparr) then
        call brdcst(msgint,n_symb,4,master,MPI_COMM_WORLD)
        call brdcst(msgint,n_k,4,master,MPI_COMM_WORLD)
        call brdcstc(s,len(s(1))*n_symb,master,MPI_COMM_WORLD)
        call brdcst(msgdbl,v_symb,8*3*n_symb,master,MPI_COMM_WORLD)
        call brdcst(msgint,i_symb,4*(n_k+1),master,MPI_COMM_WORLD)
      endif
#endif
      end subroutine read_kpath
