      subroutine boson_unpack_tau_par(f_cmpl,n1,j,f_real,n0,n)
c   ------- Unpacking from A(t) into P(t)  -------------------
c	Up  triangle of P_real =    2 * real parts of P
c	Low triangle of P_real =    2 * imaginary parts of P
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer, intent(in) :: n,n0,n1,j
	real*8, intent(in) :: f_real(n0,n0)
	complex*16, intent(out) :: f_cmpl(n1)
      integer :: i
	do i=1,j-1
	  f_cmpl(i)=0.5d0*dcmplx(f_real(i,j),f_real(j,i))
	enddo
	f_cmpl(j)=0.5d0*f_real(j,j)
	do i=j+1,n
	  f_cmpl(i)=0.5d0*dcmplx(f_real(j,i),-f_real(i,j))
	enddo
      end