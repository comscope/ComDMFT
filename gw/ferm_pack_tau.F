      subroutine ferm_pack_tau(f_cmpl,f_real,n,n0,n1)
c   ------- Packing from F(t) and F(b-t) into A and B  -------------------
c	A =     F(t) - F(b-t)
c	B =     F(t) + F(b-t)
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer, intent(in) :: n,n0,n1
	complex*16, intent(in) :: f_cmpl(n0,n0,2)
	real*8, intent(out) :: f_real(n1,n1,2)
      integer :: i,j
	do j=1,n
	  do i=1,j-1
	    f_real(i,j,1)=f_real(i,j,1)+real(f_cmpl(i,j,1))
     &		                         -real(f_cmpl(i,j,2))
	    f_real(i,j,2)=f_real(i,j,2)+real(f_cmpl(i,j,1))
     &		                         +real(f_cmpl(i,j,2))
	    f_real(j,i,1)=f_real(j,i,1)+imag(f_cmpl(i,j,1))
     &		                         -imag(f_cmpl(i,j,2))
	    f_real(j,i,2)=f_real(j,i,2)+imag(f_cmpl(i,j,1))
     &		                         +imag(f_cmpl(i,j,2))
	  enddo
	  f_real(j,j,1)=f_real(j,j,1)+real(f_cmpl(j,j,1))
     &	                           -real(f_cmpl(j,j,2))
	  f_real(j,j,2)=f_real(j,j,2)+real(f_cmpl(j,j,1))
     &	                           +real(f_cmpl(j,j,2))
	enddo
      end