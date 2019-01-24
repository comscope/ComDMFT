      subroutine boson_pack_tau(alf,f_cmpl,n0,bet,f_dble,n1,n)

c     ------- Packing from P(t) into A(tau)  -------------------
c     Up  triangle of A =   2 * dble parts of P
c     Low triangle of A =   2 * imaginary parts of P

      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: n,n0,n1
      double precision, intent(in) :: alf,bet
      complex*16, intent(in) :: f_cmpl(n0,n0)
      double precision, intent(inout) :: f_dble(n1,n1)
      integer :: i,j
      double precision :: a
      a=2.d0*alf
      if(bet==0.d0) f_dble=0.d0
      do j=1,n
        do i=1,j-1
          f_dble(i,j)=bet*f_dble(i,j)+a*dble(f_cmpl(i,j))
          f_dble(j,i)=bet*f_dble(j,i)+a*dimag(f_cmpl(i,j))
        enddo
        f_dble(j,j)=bet*f_dble(j,j)+a*dble(f_cmpl(j,j))
      enddo
      end
