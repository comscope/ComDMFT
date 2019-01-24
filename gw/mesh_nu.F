      subroutine mesh_nu
c ------- MESH for Nu ------------------------------------------------
	use manager_mod
	use parallel_mod
	use units_mod
      implicit none
	integer :: i_nu,nn,k,n
	real*8 :: pi_bet,h1,w0,wm,x,fc,fac
	pi_bet=pi/betta_t
	h1=pi_bet+pi_bet
	allocate(w_nu(0:n_nu))
	allocate(mats_num_nu(0:n_nu))
	do i_nu=0,n_nu_exa
	  w_nu(i_nu)=i_nu*h1
	enddo
c ---------- Find asymptotic part of mesh ----------------------
c ---- The part for Chebyshev's approximation ---------------------
      if((nu_max-w_nu(n_nu_exa))/n_nu_asy<10*h1) nu_max=w_nu(n_nu_exa)
     &                                                  +10*n_nu_asy*h1
      n_nu_maz=nu_max/h1
      nu_max=n_nu_maz*h1
      w_nu(n_nu)=nu_max
	allocate(cheb_nu(0:n_nu_asy-2,n_nu_asy-1))
      nn=n_nu_asy-1
      w0=sqrt(w_nu(n_nu_exa))
      wm=sqrt(nu_max)
      do i_nu=n_nu_exa+1,n_nu-1
        k=i_nu-n_nu_exa
        x=cos((0.5d0+nn-k)*pi/nn)   !! zeroes of Chebyshev's polyn.
        w_nu(i_nu)=(2.d0*w0*wm/(w0+wm-(wm-w0)*x))**2
      enddo
c ------ Calculation of CHEB_NU -------------------------------------
      do k=1,n_nu_asy-1
	  fc=(n_nu_asy-k-0.5d0)*pi/(n_nu_asy-1)
	  do n=0,n_nu_asy-2
	    fac=2.d0
	    if(n.eq.0) fac=1.d0
	    cheb_nu(n,k)=fac*cos(fc*dfloat(n))/dfloat(n_nu_asy-1)
	  enddo
	enddo	
c --- Here for compatibility when TAU--->NU -------------------------
      call mesh_nu_d			
c ------- Correspondence to true Matsubara's frequencies ---------- 
     	mats_num_nu(0)=0
	do i_nu=1,n_nu
	  mats_num_nu(i_nu)=nint(w_nu(i_nu)/h1)
	enddo
	n_nu_maz=mats_num_nu(n_nu)
	if(maswrk)write(iun,*)
     &	 ' The number of Nu-points is ',n_nu+1
	if(maswrk)write(iun,*)
     &	 ' The number of extended Nu_D-points is ',n_nu_d+1
	if(maswrk)write(iun,*)
     &	 ' The number of Matsubara Nu-frequencies is ',n_nu_maz
c ---- Now - matrices for the transform NU ----> TAU -------------
      allocate(tau_from_nu(0:n_nu,ndim3_tau,2))
	call asympt_sums_nu
	call transf_tau_from_nu_cheb
c ---------- For Interpolation purposes -------------------------      
      call for_nu_interpolation
      if(l_double_freq) then
c -------- For vertex correction to the Self Energy -------------
        allocate(tau_from_nu_d(ndim3_nu_d(me3_tau+1),0:n_tau,2))
        call transf_tau_from_nu_d
      endif
      call timel('**** MESH_NU finished **************')
      end
