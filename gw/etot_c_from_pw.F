      subroutine etot_c_from_pw(n_pbt,iq,i_nu,p,ws)
      use atom_mod
      use etot_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
      integer, intent(in) :: n_pbt,iq,i_nu
	complex*16, intent(in) :: p(n_pbt,n_pbt),ws(n_pbt,n_pbt)
	integer :: j,i
	real*8 :: de
	de=0.d0
	do j=1,n_pbt
	  do i=1,n_pbt
	    de=de-p(i,j)*ws(j,i)
	  enddo
	enddo
      e_c_solid_pw=e_c_solid_pw+et_sum_nu(i_nu)*wgt(iq)*de
      end