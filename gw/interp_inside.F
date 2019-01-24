      subroutine interp_inside(f_in,f_out,k_in,k_out,n_in,n_out)
c ------- Lagrange's interpolation ----------
c     The boundary points coincide
	use manager_mod
	use parallel_mod
	use solid_mod
      implicit none
      integer, intent(in) :: n_in,n_out
	real*8, intent(in) :: f_in(0:n_in),k_in(0:n_in),k_out(0:n_out)
	real*8, intent(out) :: f_out(0:n_out)
	integer :: i,j,k,nprev,i0
	nprev=0
	do i=0,n_out
c --- Renew NPREV -------------------
        k=nprev
        do j=nprev+1,n_in-1
          if(k_in(j)>k_out(i)) exit
          k=j
        enddo
        nprev=k
        i0=nprev-2
        i0=max(i0,0)
        i0=min(i0,n_in-5)
        call fint(f_in(i0),k_in(i0),6,k_out(i),f_out(i))
      enddo
      end