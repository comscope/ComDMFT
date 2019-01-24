      subroutine g_x_omega(ispin,k,gx_omega,w_n,ndim,x)
c ------ Exchange part of Green's Function is calculated here -------------
c	x - Chemical Potential
	use manager_mod
	use solid_mod
      implicit none
	integer, intent(in) :: ndim
	real*8, intent(in) :: w_n
	complex*16, intent(out) :: gx_omega(ndim)
	integer :: i,k,ispin
	real*8 :: de,x
	complex*16 :: cf
	do i=1,ndim
        de=e_bnd(i,k,ispin)-x
	  cf=dcmplx(-de,w_n)
	  cf=(1.d0,0.d0)/cf
	  gx_omega(i)=cf
	enddo
      end