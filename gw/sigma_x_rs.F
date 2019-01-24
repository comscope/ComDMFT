      subroutine sigma_x_rs(ispin,v_rs_mm)
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer, intent(in) :: ispin
      real*8, intent(in) :: v_rs_mm(n_pbmtm,n_pbmtm,ndim_nrr(me3_k+1))
	real*8 :: dum
	complex*16 :: dumc
	real*8, allocatable :: g_r(:,:,:)
	complex*16, allocatable :: s_c(:,:,:),s_c1(:,:,:,:),g_c1(:,:,:,:)
c
c	---------- MT_MT Part ------------------------------------
c
	allocate(s_c(nfun,nfun,ndim3_kk(me3_k+1)))
	call sigx_mm_r(ispin,s_c,v_rs_mm)
	call sigx_mm_k_from_r(s_c)
	deallocate(s_c)
      call timel('* Sigma_x_rs: MM finished **********')
c
c	---------- INT_MT Part ------------------------------------
c
	if(irel<=1) then
	  allocate(g_r(nr_red,nfun,ndim3_kk(me3_k+1)))
	  call gx_rs_from_ks_mi(g_r,dumc,ispin)
	else if(irel==2) then
	  allocate(g_c1(nr_red,nfun,ndim3_kk(me3_k+1),nrel))
	  call gx_rs_from_ks_mi(dum,g_c1,ispin)
	endif
	allocate(s_c1(nr_red,nfun,nrel,ndim3_kk(me3_k+1)))
	if(irel<=1) then
	  call sigx_im_r(ispin,g_r,s_c1)
	  deallocate(g_r)
	else if(irel==2) then
	  call sigx_im_r_rel(ispin,g_c1,s_c1)
	  deallocate(g_c1)
	endif
	call sigx_im_k_from_r(s_c1)
	deallocate(s_c1)
      call timel('* Sigma_x_rs: MI finished **********')
c
c	---------- INT_INT Part ------------------------------------
c
	call sigx_ii_k_from_r(ispin)
      call timel('* Sigma_x_rs: II finished **********')
      end