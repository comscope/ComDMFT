      subroutine from_tau_r_to_tau_herm(f_tau_h,f_tau_r,n0,n1,n)
c   ------- Transform TAU ---> NU representation --------------------------
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer, intent(in) :: n,n0,n1
	real*8, intent(in) :: f_tau_r(n1,n1,2,ndim3_tau)
	real*8, intent(out) :: f_tau_h(n0,n0,2,ndim3_tau)
      integer :: i,j,ind_tau
	f_tau_h=0.d0
	do ind_tau=1,ndim3_tau
	  do j=1,n
	    do i=1,j
	      f_tau_h(i,j,1,ind_tau)=f_tau_r(i,j,1,ind_tau)
     &	                          -f_tau_r(i,j,2,ind_tau)
	      f_tau_h(i,j,2,ind_tau)=f_tau_r(i,j,1,ind_tau)
     &	                          +f_tau_r(i,j,2,ind_tau)
	    enddo
	  enddo
	enddo
      end