      subroutine invers_h(n,a,n0)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: n,n0
      complex*16, intent(inout) :: a(n0,n0)
      integer :: info,j,i
      integer, allocatable :: ipiv(:)
      complex*16,allocatable :: wrk(:,:)
      allocate(ipiv(n))
      allocate(wrk(n,n))
      call zhetrf('U',n,a,n0,ipiv,wrk,n*n,info)
      call zhetri('U',n,a,n0,ipiv,wrk,info)
      do j=1,n
        do i=j+1,n
          a(i,j)=conjg(a(j,i))
        enddo
      enddo
      deallocate(ipiv,wrk)
      end
