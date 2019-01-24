      subroutine w_red_mi(wmi,w_v,npl)
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
      implicit none
	integer, intent(in) :: npl
	real*8, intent(in) :: w_v(n_pbtot,n_pbtot)
	complex*16, intent(out) :: wmi(n_pbmt,npl)
	integer :: n0
	complex*16, allocatable :: tmp(:,:)
      n0=n_pbmt+npl
	allocate(tmp(n0,n0))
	call boson_unpack_tau(tmp,n0,w_v,n_pbtot,n0)
	wmi=tmp(1:n_pbmt,n_pbmt+1:n_pbmt+npl)
	deallocate(tmp)
      end