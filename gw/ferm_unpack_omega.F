      subroutine ferm_unpack_omega(f_cmpl,f_real,n,n0,n1)
c   ------- Unpacking from A(w) and B(w) into F(w)  -------------------
c	F(w) = A(w) + i B(w)
c     But A and B are stored in square R*8 form
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer, intent(in) :: n,n0,n1
	real*8, intent(in) :: f_real(n0,n0,2)
	complex*16, intent(out) :: f_cmpl(n1,n1)
      integer :: i,j
	do j=1,n
	  do i=1,j-1
	    f_cmpl(i,j)=dcmplx(f_real(i,j,1)-f_real(j,i,2),
     &		                 f_real(j,i,1)+f_real(i,j,2))
	  enddo
	  f_cmpl(j,j)=dcmplx(f_real(j,j,1),f_real(j,j,2))
	  do i=j+1,n
	    f_cmpl(i,j)=dcmplx(f_real(j,i,1)+f_real(i,j,2),
     &		                 f_real(j,i,2)-f_real(i,j,1))
	  enddo
	enddo
      end