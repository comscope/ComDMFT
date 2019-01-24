      subroutine mesh_tau_cheb
	use manager_mod
	use parallel_mod
	use units_mod
      implicit none
	integer :: k,i_tau,l
	real*8 :: pi_n,fc,fac,xk
	tau_mesh(0)=0.d0
	tau_mesh(n_tau)=betta_t
c ---------------------------------------------------------------------
      pi_n=pi/dfloat(2*(n_tau-1))
	do i_tau=1,n_tau/2
	  xk=cos(pi_n*dfloat(2*(n_tau-i_tau)-1))  !! Cheb. roots
	  tau_mesh(i_tau)=0.5d0*betta_t*(1.d0+xk)**3
	  tau_mesh(n_tau-i_tau)=betta_t-tau_mesh(i_tau)
	enddo
	allocate(cheb_tau_even(0:n_tau/2-1,n_tau/2),
     &         cheb_tau_odd(0:n_tau/2-2,n_tau/2-1))
c ------ Calculation of CHEB_TAU_EVEN ---------------------------------
      do k=1,n_tau/2-1
	  fc=dfloat(n_tau-k-k)/(n_tau-1)*pi
	  do l=0,n_tau/2-1
	    fac=2.d0
	    if(l==0) fac=1.d0
	    cheb_tau_even(l,k)=2.d0*(-1)**l*fac/(n_tau-1)*cos(fc*l)
	  enddo
	enddo
	do l=0,n_tau/2-1
	  fac=2.d0
	  if(l==0) fac=1.d0
	  cheb_tau_even(l,n_tau/2)=(-1)**l*fac/(n_tau-1)
	enddo
c ------ Calculation of CHEB_TAU_ODD ----------------------------------
      do k=1,n_tau/2-1
	  fc=dfloat(n_tau-k-k)/2.d0/(n_tau-1)*pi
	  fac=4.d0/(n_tau-1)
	  do l=0,n_tau/2-2
	    cheb_tau_odd(l,k)=fac*(-1)**(l+1)*sin(fc*(l+l+1))
	  enddo
	enddo
      end