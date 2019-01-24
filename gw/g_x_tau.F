      subroutine g_x_tau(ispin,k,gx_tau,t,ndim,x)
c ------ Exchange part of Green's Function is calculated here -------------
	use manager_mod
	use solid_mod
      implicit none
	integer, intent(in) :: ispin,k,ndim
	real*8, intent(in) :: t,x
	real*8, intent(out) :: gx_tau(ndim)
	integer :: i
	real*8 :: de,green0
c ------------------------------------------------------------
	do i=1,ndim
	  de=e_bnd(i,k,ispin)-x
	  gx_tau(i)=green0(de,t)
	enddo
      end