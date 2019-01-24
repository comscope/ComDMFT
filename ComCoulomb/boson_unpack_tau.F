      subroutine boson_unpack_tau(f_cmpl,n1,f_dble,n0,n)
c     ------- Unpacking from A(t) into P(t)  -------------------
c     Up  triangle of P_dble =    2 * dble parts of P
c     Low triangle of P_dble =    2 * imaginary parts of P
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: n,n0,n1
      double precision, intent(in) :: f_dble(n0,n0)
      complex*16, intent(out) :: f_cmpl(n1,n1)
      integer :: i,j
      do j=1,n
        do i=1,j-1
          f_cmpl(i,j)=0.5d0*dcmplx(f_dble(i,j),f_dble(j,i))
          f_cmpl(j,i)=conjg(f_cmpl(i,j))
        enddo
        f_cmpl(j,j)=0.5d0*f_dble(j,j)
      enddo
      end
