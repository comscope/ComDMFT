      subroutine ferm_pack_omega(f_cmpl,f_real,n,n0,n1)
c   ------- Packing from C*16 F(w) into A(w) and B(w)  -------------------
c	A(w) =  1/2 (F(w) + F^(+)(w))  We store them in R*8 square forms
c	B(w) = -i/2 (F(w) - F^(+)(w))
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer, intent(in) :: n,n0,n1
	complex*16, intent(in) :: f_cmpl(n0,n0)
	real*8, intent(out) :: f_real(n1,n1,2)
      integer :: i,j
	real*8 :: a,b,c,d
	do j=1,n
	  do i=1,j-1
	    a=real(f_cmpl(i,j))/2
	    b=imag(f_cmpl(i,j))/2
	    c=real(f_cmpl(j,i))/2
	    d=imag(f_cmpl(j,i))/2
	    f_real(i,j,1)=a+c
	    f_real(i,j,2)=b+d
	    f_real(j,i,1)=b-d
	    f_real(j,i,2)=-a+c
	  enddo
	  f_real(j,j,1)=real(f_cmpl(j,j))
	  f_real(j,j,2)=imag(f_cmpl(j,j))
	enddo
      end