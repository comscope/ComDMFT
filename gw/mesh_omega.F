      subroutine mesh_omega
c ------- MESH for Omega ------------------------------------------------
	use manager_mod
	use parallel_mod
	use units_mod
      implicit none
	integer :: i_omega,n,nn,k
	real*8 :: pi_bet,h1,w0,wm,x,fc,fac
	pi_bet=pi/betta_t
	h1=pi_bet+pi_bet
	allocate(w_omega(0:n_omega))
	allocate(mats_num_omega(0:n_omega))
	do i_omega=0,n_omega_exa
	  w_omega(i_omega)=pi_bet+i_omega*h1
	enddo
c ---------- Find asymptotic part of mesh -------------------------
c ---- The part for Chebyshev's approximation ---------------------
      if((omega_max-w_omega(n_omega_exa))/n_omega_asy<10*h1) 
     &  omega_max=w_omega(n_omega_exa)+10*n_omega_asy*h1
      n_omega_maz=(omega_max-pi_bet)/h1
      omega_max=(2*n_omega_maz+1)*pi_bet
      w_omega(n_omega)=omega_max
	allocate(cheb_omega(0:n_omega_asy-2,n_omega_asy-1))
      nn=n_omega_asy-1
      w0=sqrt(w_omega(n_omega_exa))
      wm=sqrt(omega_max)
      do i_omega=n_omega_exa+1,n_omega-1
        k=i_omega-n_omega_exa
        x=cos((0.5d0+nn-k)*pi/nn)
        w_omega(i_omega)=(2.d0*w0*wm/(w0+wm-(wm-w0)*x))**2
      enddo
c ------ Calculation of CHEB_OMEGA -------------------------------------
      do k=1,n_omega_asy-1
	  fc=(n_omega_asy-k-0.5d0)*pi/(n_omega_asy-1)
	  do n=0,n_omega_asy-2
	    fac=2.d0
	    if(n.eq.0) fac=1.d0
	    cheb_omega(n,k)=fac*cos(fc*dfloat(n))/dfloat(n_omega_asy-1)
	  enddo
	enddo
c --- Here for compatibility when TAU--->OMEGA -------------------------
      call mesh_omega_d
c ------- Correspondence to true Matsubara's frequencies ---------- 
     	mats_num_omega(0)=0
	do i_omega=1,n_omega
	  mats_num_omega(i_omega)=nint((w_omega(i_omega)-pi_bet)/h1)
	enddo
	n_omega_maz=mats_num_omega(n_omega)
	if(maswrk)write(iun,*)
     &	 ' The number of Omega-points is ',n_omega+1
	if(maswrk)write(iun,*)
     &	 ' The number of extended Omega_D-points is ',n_omega_d+1
	if(maswrk)write(iun,*)
     &	 ' The number of Matsubara Omega-frequencies is ',n_omega_maz
c ---- Now - matrices for the transform OMEGA ----> TAU -------------
      allocate(tau_from_omega(0:n_omega,ndim3_tau,2))
	call asympt_sums_omega
	call transf_tau_from_omega_cheb
c ---------- For Interpolation purposes -------------------------
      call for_omega_interpolation
      call timel('**** MESH_OMEGA finished ***********')
      end