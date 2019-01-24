      subroutine pack_hermit(a_compl,a_real,n,n0,n1,alpha,beta)
      implicit none
	integer, intent(in) :: n,n0,n1
	real*8, intent(in) :: alpha,beta
	complex*16, intent(in) :: a_compl(n0,n0)
	real*8, intent(out) :: a_real(n1,n1)
	integer :: i,j
	if(alpha==0.d0) a_real=0.d0
	do j=1,n
	  do i=1,j-1
	    a_real(i,j)=alpha*a_real(i,j)+beta*real(a_compl(i,j))
	    a_real(j,i)=alpha*a_real(j,i)+beta*imag(a_compl(i,j))
	  enddo
	  a_real(j,j)=alpha*a_real(j,j)+beta*real(a_compl(j,j))
	enddo
      end
      
      
      
	subroutine unpack_hermit(a_compl,a_real,n,n0,n1,alpha,beta)
      implicit none
	integer, intent(in) :: n,n0,n1
	complex*16, intent(in) :: alpha,beta
	real*8, intent(in) :: a_real(n0,n0)
	complex*16, intent(out) :: a_compl(n1,n1)
	integer :: i,j
	if(alpha==(0.d0,0.d0)) a_compl=(0.d0,0.d0)
	do j=1,n
	  do i=1,j-1
	    a_compl(i,j)=alpha*a_compl(i,j)
     &		+beta*dcmplx(a_real(i,j),a_real(j,i))
	    a_compl(j,i)=conjg(a_compl(i,j))
	  enddo
	  a_compl(j,j)=alpha*a_compl(j,j)+beta*a_real(j,j)
	enddo
      end
      
      subroutine pack_sym(a_compl,a_real,n,n0,n1,alpha,beta)
      implicit none
	integer, intent(in) :: n,n0,n1
	real*8, intent(in) :: alpha,beta
	real*8, intent(in) :: a_compl(n0,n0)
	real*8, intent(out) :: a_real(n1,n1)
	integer :: i,j
	if(alpha==0.d0) a_real=0.d0
	do j=1,n
	  do i=1,j
	    a_real(i,j)=alpha*a_real(i,j)+beta*a_compl(i,j)
	  enddo
	enddo
      end