      subroutine sigma_prelim(ispin,v_rs_mm)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
	use units_mod
	implicit none
	integer, intent(in) :: ispin
      real*8, intent(in) :: v_rs_mm(n_pbmtm,n_pbmtm,ndim_nrr(me3_k+1))
      pv_x=(0.d0,0.d0)
      call get_pc(ispin)
      call timel('******* get_pc finished ************')
c	call get_pv(ispin)
c      call timel('******* get_pv finished ************')
      allocate(pvv(nfun,nfun,maxval(ndim3_k)))
      allocate(pvg(nfun,nbasmpw,maxval(ndim3_k)))
	pvv=(0.d0,0.d0)
	pvg=(0.d0,0.d0)
	allocate(pgg(nbasmpw,nbasmpw,maxval(ndim3_k)))
	pgg=(0.d0,0.d0)
	call sigma_x_rs(ispin,v_rs_mm)
      call timel('**** SIGMA_PRELIM finished *********')
      end
