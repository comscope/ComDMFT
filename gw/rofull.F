      SUBROUTINE ROfull(correlation,key)
**********************************************************
* The calculation of the  valence density                *
**********************************************************
c     key = 0  - No density normalization control
c     key = 1  - We renormalize the density
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
      IMPLICIT none
      logical, intent(in) :: correlation
      integer, intent(in) :: key
      logical :: jreal,jimag
	integer :: keym

      rointr=(0.d0,0.d0)
      if(magn.eq.2) spintr=(0.d0,0.d0)


! ro_int_g calculates the interstitial sector of the charge density, rointr, and symmetrizes it.
! If doing irel<=1, 
! rointr = \sum_k ev_bnd(k) * g * ev_bnd^\dagger * w(k) / amega, plus symmetrization.
! w(k) = wgt(k)*2.d0/dfloat(nrel*nspin)
! ev_bnd is the interstitial sector of the DFT bands.
! g is g_full_00 if doing dft/hf but not with dmft, and is g_full_0 if doing dmft or qp or gw.
! If doing irel=2, the equations are a bit different.
	call ro_int_g(correlation)

	  ro=0.d0
        if(magn.eq.2) spmt=0.d0

! The get_t+ro_val pair calculates the muffin-tin sector of the charge density.
        call get_t
	  call ro_val

	if(key==1) call renr_val
	  call add_cor(1)
        if(itest.gt.0) then
          keym=0
          if(nspin.eq.2) keym=2
	    text=' Matching : New Density!  '
          call matching(keym,ro,rointr)
          if(nspin.eq.2) then
	      text=' Matching : New Density of Spin Moment - Z!'
            call matching(3,ro,rointr)
	    endif
	    if(magn.eq.2) then
            jreal=.false.
            jimag=.false.
            if(inv.lt.2) jreal=.true.
            if(inv.ne.1) jimag=.true.
	      text=' Matching : New Density of Spin Moment!'
            call matchav(spmt,spintr,jreal,jimag)
          endif
        endif  !!! over itest
	  call mag_mom
      call timel('****** Rofull finished *************')
      END
