      subroutine sigma_x_int_int(nbas,ndim,sig_x0,pii)
	use atom_mod
	use manager_mod
	use solid_mod
      implicit none
      integer, intent(in) :: nbas,ndim
	complex*16, intent(in) :: pii(nbasmpw,nbasmpw)
      complex*16 :: sig_x0(ndim,ndim)
	sig_x0(nlcloc+1:ndim,nlcloc+1:ndim)=
     &	sig_x0(nlcloc+1:ndim,nlcloc+1:ndim)+pii(1:nbas,1:nbas)
      end