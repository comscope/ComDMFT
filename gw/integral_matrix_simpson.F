      subroutine integral_matrix_simpson(a,b,x,k,m0,n0,m1,n1,m,n,alpha)
      implicit none
	integer, intent(in) :: k,m0,n0,m1,n1,m,n
      real*8, intent(in) :: a(m0,n0,0:k),x(0:k),alpha
      real*8, intent(out) :: b(m1,n1)
	integer :: i,j
	real*8, allocatable :: w(:)
	if(k==0) then
	  b=alpha*b
	  return
	endif
	allocate(w(0:k))
	call simpson_ihmg_weights(x,k,w)
	do j=1,n
	  do i=1,m
	    b(i,j)=alpha*b(i,j)+dot_product(a(i,j,0:k),w(0:k))
	  enddo
	enddo
	deallocate(w)
      end