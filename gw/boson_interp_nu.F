      subroutine boson_interp_nu(n,a_in,a_cheb,nu_mats,a_out)
c ---- Interpolation A (bosonic) onto arbitrary Matsubara' NU ----
c --- Input and output arrays are in compact form (triangle if needed)
c --- where n - total number of elements ------------------------------
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
      integer, intent(in) :: n,nu_mats
	real*8, intent(in) :: a_in(n,0:n_nu),
     &                      a_cheb(n,0:n_nu_asy-2)
	real*8, intent(out) :: a_out(n)
	integer :: i,nua_mats
	real*8 :: x,x_cheb_from_nu,v2,v4,pi8,bt,v
	real*8, allocatable :: cheb(:)
	pi8=8.d0*pi
	bt=(pi+pi)/betta_t
	nua_mats=abs(nu_mats)
	v=nua_mats*bt
      if(nua_mats<=n_nu_exa) then
        a_out=a_in(:,nua_mats)
      else if(nua_mats<=n_nu_maz) then
	  allocate(cheb(0:n_nu_asy-2))
        x=x_cheb_from_nu(v)
	  call dcheb_recurr(x,cheb,n_nu_asy-2)
	  a_out=0.d0
	  do i=0,n_nu_asy-2
	    a_out=a_out+a_cheb(:,i)*cheb(i)
	  enddo
	  deallocate(cheb)
      else
        v2=v*v
        v4=v2*v2
	  a_out=(asy_nu_interp(1,1,1)*a_in(:,n_nu-1)
     &	      +asy_nu_interp(2,1,1)*a_in(:,n_nu))/v2
     &       +(asy_nu_interp(1,2,1)*a_in(:,n_nu-1)
     &	      +asy_nu_interp(2,2,1)*a_in(:,n_nu))/v4
	endif
      end