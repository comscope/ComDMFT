      subroutine print_matrix(iunit,x,ndx,m,n)
c
c Print a matrix to standard output.
c
c The matrix is printed in nmax columns at a time.
c Arguments:
c iunit - the output unit number
c x     - the matrix in double precision
c ndx   - the leading dimension of x
c m     - the number of rows
c n     - the number of columns
c
      implicit none
      integer iunit
      integer ndx
      integer m, n
      real*8 x(ndx,n)
      integer i,j,nm,nmax
      parameter(nmax=8)
      do nm = 1, n, nmax
        write(iunit,'(4x,8i12)')(j,j=nm,min(n,nm+nmax-1))
        do i = 1, m
          write(iunit,'(i3,x,8f12.6)')i,(x(i,j),j=nm,min(n,nm+nmax-1))
        enddo
        write(iunit,*)
      enddo
      end
