      subroutine read_kpath(mx_symb,mx_k,n_symb,num_k_path,s,
     $  v_symb,ind_symb)
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
c     X   0.00  -0.50   0.00       G   0.00   0.00   0.00
c     Y   0.50   0.00   0.00       L   0.50  -0.50   0.00
c     G   0.00   0.00   0.00       Z  -0.50   0.00   0.50
c     N  -0.50  -0.50   0.50       G   0.00   0.00   0.00
c     M   0.00   0.00   0.50       R   0.00  -0.50   0.50
c     G   0.00   0.00   0.00
c     
      integer, intent(in) ::  mx_symb ! the max number of high symmetry
! points
      integer, intent(in) ::  mx_k ! the max number of k-points on
! the path
      integer, intent(out) :: n_symb ! the number of high symmetry 
!     points
      integer, intent(out) :: num_k_path ! the number of k-points on the
! path
      character*3, intent(out) :: s(mx_symb) ! the high symmetry
! point labels
      double precision, intent(out) :: v_symb(3,mx_symb) ! the k-point
! coordinates
      integer, intent(out) :: ind_symb(mx_k) ! the k-point indexes on
! the path
c     
      integer :: sym_check,isymb,nlines,ii,ik
      integer :: ierr,io           ! error code
      double precision :: tol = 2.0d-8 ! tolerance
      character*3      :: s1    ! high symmetry point label as read
      double precision :: v1(3) ! k-point coordinates as read      
c     

      
      if (maswrk) then

c$$$        write(iun,'(a)') 'read_kpath, here'
        nlines = 0 
        OPEN (1, file = 'kpath')
        DO
          READ(1,*,iostat=io)
          IF (io/=0) EXIT
          nlines = nlines + 1
        END DO
        CLOSE (1)      
        
        open(3,file="kpath",status="old",form="formatted")

        num_k_path=0
        n_symb=0
        s        = ' '
        v_symb   = 0.0d0
        ind_symb = 0        

        do ii=1, nlines
          call inp_line(3,ierr)
          if (ierr.eq.ERR_EOF) exit          
          do ik=1, 2
            call inp_str(s1,ierr)
            call inp_dbl(v1(1),ierr)
            call inp_dbl(v1(2),ierr)
            call inp_dbl(v1(3),ierr)
            
c$$$            write(iun,*) s1, v1
            
            if (
     $        (num_k_path .eq. 0) .or.
     $        ((num_k_path .gt. 0) .and.
     $        (s(ind_symb(num_k_path)) .ne. s1))
     $        ) then

c$$$              write(iun,'(i5, 100a3)') num_k_path, s(1:num_k_path)
c$$$              write(iun,*) s1, s(ind_symb(num_k_path))              

              num_k_path=num_k_path+1
              
              sym_check=0
              do isymb = 1, n_symb
                if (s(isymb).eq.s1.and.
     &            abs(v_symb(1,isymb)-v1(1)).lt.tol.and.
     &            abs(v_symb(2,isymb)-v1(2)).lt.tol.and.
     &            abs(v_symb(3,isymb)-v1(3)).lt.tol) then
                  ind_symb(num_k_path)=isymb
                  sym_check=sym_check+1
                endif
              enddo
              
              if (sym_check.eq. 0) then
                n_symb=n_symb+1
                s(n_symb)=s1
                v_symb(:,n_symb)=v1
                ind_symb(num_k_path)=n_symb 
              endif
            endif
          enddo
        enddo
        
        
c$$$        write(iun,'(a, 2i5)') 'kpath', num_k_path, n_symb
c$$$        do ii=1, n_symb
c$$$          write(iun,'(a, 3f12.6)') s(ii), v_symb(:,ii)
c$$$        enddo
c$$$        write(iun,'(40i5)') (ind_symb(ii), ii=1, num_k_path)
      endif
#ifdef MPI
      if (goparr) then
        call brdcst(msgint,n_symb,4,master,MPI_COMM_WORLD)
        call brdcst(msgint,num_k_path,4,master,MPI_COMM_WORLD)
        call brdcstc(s,len(s(1))*mx_symb,master,MPI_COMM_WORLD)
        call brdcst(msgdbl,v_symb,8*3*mx_symb,master,MPI_COMM_WORLD)
        call brdcst(msgint,ind_symb,4*mx_k,
     $    master,MPI_COMM_WORLD)
      endif
#endif
      end subroutine read_kpath
      
