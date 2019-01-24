      real*8 function check_hermit(a,n,n0)
	implicit none
	integer :: n,i,j,n0
	real*8 :: d
	complex*16 :: a(n0,n0)
	d=0.d0
	do j=1,n
	  do i=1,j
	    d=max(d,abs(a(i,j)-conjg(a(j,i))))
	    a(i,j)=0.5d0*(a(i,j)+conjg(a(j,i)))
	    a(j,i)=conjg(a(i,j))
	  enddo
	enddo
	check_hermit=d
      end
      
      subroutine check_sym_matr(a,n,n0,err)
	implicit none
	integer, intent(in) :: n,n0
	real*8, intent(inout) :: a(n0,n0)
	real*8, intent(out) :: err
	integer :: i,j
	err=0.d0
	do j=1,n
	  do i=1,j-1
	    err=max(err,abs(a(i,j)-a(j,i)))
	    a(i,j)=0.5d0*(a(i,j)+a(j,i))
	    a(j,i)=a(i,j)
	  enddo
	enddo
      end