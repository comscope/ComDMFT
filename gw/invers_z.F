	subroutine invers_z(n,a,n0)
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer, intent(in) :: n,n0
  	complex*16, intent(inout) :: a(n0,n0)
	integer :: info
	integer, allocatable :: ipiv(:)
	complex*16,allocatable :: wrk(:,:)
	allocate(ipiv(n))
	allocate(wrk(n,n))
	call zgetrf(n,n,a,n0,ipiv,info)
	call zgetri(n,a,n0,ipiv,wrk,n*n,info)
	deallocate(ipiv,wrk)
      end
      
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
      
      
      subroutine invers_r(n,a,n0)
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer, intent(in) :: n,n0
  	real*8, intent(inout) :: a(n0,n0)
	integer :: info
	integer, allocatable :: ipiv(:)
	real*8,allocatable :: wrk(:,:)
	allocate(ipiv(n))
	allocate(wrk(n,n))
	call dgetrf(n,n,a,n0,ipiv,info)
	call dgetri(n,a,n0,ipiv,wrk,n*n,info)
	deallocate(ipiv,wrk)
      end
      
      subroutine sqrt_matr(n,n0,n1,a,b)
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
      implicit none
      integer, intent(in) :: n,n0,n1
      real*8, intent(in) :: a(n0,n0)
      real*8, intent(out) :: b(n1,n1)
      integer :: j,i,il
	real*8, allocatable :: c(:,:),e(:)
	allocate(c(n,n))
	allocate(e(n))
	c=a(1:n,1:n)
      call eig_val_solver_r(n,n,c,e)
	b=0.d0
	do j=1,n
	  do i=1,n
	    do il=1,n
	      b(i,j)=b(i,j)+c(i,il)*sqrt(e(il))*c(j,il)
	    enddo
	  enddo
	enddo
	deallocate(c,e)
      end