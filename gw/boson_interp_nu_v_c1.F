      subroutine boson_interp_nu_v_c1(n,a_in,a_cheb,nu,a_out)
c ---- Interpolation A (bosonic) onto arbitrary NU ----
c --- Input and output arrays are in compact form (triangle if needed)
c --- where n - total number of elements ------------------------------
c -------- Assume f(-v)=f*(v) -------------------------------------
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
      integer, intent(in) :: n
	real*8, intent(in) :: nu
	complex*16, intent(in) :: a_in(n,0:n_nu),a_cheb(n,0:n_nu_asy-2)
	complex*16, intent(out) :: a_out(n)
	integer :: i,nua_mats
	real*8 :: x,x_cheb_from_nu,v2,v4,v,pi2b
	real*8, allocatable :: cheb(:)
	pi2b=(pi+pi)/betta_t
	v=abs(nu)
      if(v<=w_nu(n_nu_exa)+1.d-12) then
        nua_mats=idnint(v/pi2b)
        a_out=a_in(:,nua_mats)
      else if(v<=w_nu(n_nu)+1.d-12) then
	  allocate(cheb(0:n_nu_asy-2))
        x=x_cheb_from_nu(v)
	  call dcheb_recurr(x,cheb,n_nu_asy-2)
	  a_out=(0.d0,0.d0)
	  do i=0,n_nu_asy-2
	    a_out=a_out+a_cheb(:,i)*cheb(i)
	  enddo
	  deallocate(cheb)
      else
        v2=w_nu(n_nu)/v
        v4=v2*v2
        do i=1,n
          a_out(i)=dcmplx(real(a_in(i,n_nu))*v4,
     &                    imag(a_in(i,n_nu))*v2)
        enddo
	endif
	if(nu<0.d0) a_out=conjg(a_out)
      end