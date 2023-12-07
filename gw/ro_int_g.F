! ro_int_g calculates the interstitial sector of the charge density, rointr, and symmetrizes it.
! If doing irel<=1, 
! rointr = \sum_k ev_bnd(k) * g * ev_bnd^\dagger * w(k) / amega, plus symmetrization.
! w(k) = wgt(k)*2.d0/dfloat(nrel*nspin)
! ev_bnd is the interstitial sector of the DFT bands.
! g is g_full_00 if doing dft/hf but not with dmft, and is g_full_0 if doing dmft or qp or gw.
! If doing irel=2, the equations are a bit different.
      SUBROUTINE ro_int_g(correlation)
	use manager_mod
	use parallel_mod
	use solid_mod
      IMPLICIT none
      logical, intent(in) :: correlation
	integer :: ispin,k,ind_k
	logical :: jreal,jimag
      DO ISPIN=1,NSPIN
	  do ind_k=1,ndim3_k(me3_k+1)
	    k=n3_mpi_k(me3_k+1)+ind_k

! sumi_gw calculates the interstitial sector of the charge density, rointr.
! rointr is later symmetrized by ro_int_g.
! If doing irel<=1, 
! rointr = \sum_k ev_bnd(k) * g * ev_bnd^\dagger * w(k) / amega
! w(k) = wgt(k)*2.d0/dfloat(nrel*nspin)
! ev_bnd is the interstitial sector of the DFT bands.
! g is g_full_00 if doing dft/hf but not with dmft, and is g_full_0 if doing dmft or qp or gw.
! If doing irel=2, the equations are a bit different.
	    call sumi_gw(correlation,ind_k,k,ispin)
        enddo 
        if(nproc_k/=1) then
	    call DGOP(rointr(1,ispin),2*nplwro,'  +',comm_pnt)
          if(magn.eq.2) then
	      call DGOP(spintr,6*nplwro,'  +',comm_pnt)
	    endif
	  endif
c ------- Symmetrization ----------------------------------------------
        call symscal(rointr(1,ispin))
      enddo   !!! over ispin
      if(magn.eq.2) then
        jreal=.false.
        jimag=.false.
        if(inv.lt.2) jreal=.true.
        if(inv.ne.1) jimag=.true.
        call symvec(spintr,jreal,jimag)
      endif  !!! for magn = 2 only
      END
