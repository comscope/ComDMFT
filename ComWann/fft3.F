      
      subroutine fft3(n1,n2,n3,nx,f,key)
c     key=1  - from K to R
c     key=-1 - from R to K
      use comwann_mod
      implicit none
      integer, intent(in) :: n1,n2,n3,nx,key
      double precision, intent(inout) :: f(2*n1*n2*n3,nx)
      integer :: j
c     --------------------------------------------------------
      do j=1,nx
        call cft(f(1,j),f(2,j),n1*n2*n3,n1,n1,2*key)
        call cft(f(1,j),f(2,j),n1*n2*n3,n2,n1*n2,2*key)
        call cft(f(1,j),f(2,j),n1*n2*n3,n3,n1*n2*n3,2*key)
      enddo
      end
