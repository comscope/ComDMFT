      subroutine boson_pack_tau(alf,f_cmpl,n0,bet,f_real,n1,n)
c   ------- Packing from P(t) into A(tau)  -------------------
c	Up  triangle of A =   2 * real parts of P
c	Low triangle of A =   2 * imaginary parts of P
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer, intent(in) :: n,n0,n1
	real*8, intent(in) :: alf,bet
	complex*16, intent(in) :: f_cmpl(n0,n0)
	real*8, intent(inout) :: f_real(n1,n1)
      integer :: i,j
      real*8 :: a
      a=2.d0*alf
      if(bet==0.d0) f_real=0.d0
	do j=1,n
	  do i=1,j-1
	    f_real(i,j)=bet*f_real(i,j)+a*real(f_cmpl(i,j))
	    f_real(j,i)=bet*f_real(j,i)+a*imag(f_cmpl(i,j))
	  enddo
	  f_real(j,j)=bet*f_real(j,j)+a*real(f_cmpl(j,j))
	enddo
      end