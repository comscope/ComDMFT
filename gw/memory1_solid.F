      subroutine memory1_solid
	use atom_mod
	use etot_mod
	use manager_mod
	use parallel_mod
	use solid_mod
      implicit none
	integer :: me1,nlmb
	me1=me+1
	nlmb=nrel*maxb+1
	allocate(g_loc_0(maxel,maxel,natom,nspin_0))
      allocate(ro(maxmt))
      allocate(p_f(0:maxnrad,maxntle,nlmb,nsort,nspin))
      allocate(q_f(0:maxnrad,maxntle,nlmb,nsort,nspin))
      allocate(gfun_old(maxwf,nspin))
	allocate(gfund_old(maxwf,nspin))
	allocate(t_ro(maxlfun,maxlfun,maxnsym,nsort,nspin_0))
      allocate(v_mt_h(maxmt/nspin_0))
      allocate(v_mt_h_old(maxmt/nspin_0))
      allocate(v_mt_xc(maxmt))
      v_mt_xc=0.d0
      if(rad_solver/='DFT') then
        allocate(f_mt_ex(0:maxnrad,maxntle,nlmb,nsort,nspin))
        allocate(f_mt_exq(0:maxnrad,maxntle,nlmb,nsort,nspin))
      endif
      allocate(pv(maxel,maxel,natom,nspin_0))
      allocate(pvj(maxel,maxel,natom))
	allocate(pv_x(maxel,maxel,natom))
	allocate(v_mm(nlb2,nlb2,maxval(ndim3_k)))
	allocate(g_full_00(nbndf,maxval(ndim3_k),nspin_0))
      if(iter_gw+iter_qp/=0) then
	  allocate(g_full_0(nbndf,nbndf,maxval(ndim3_k),nspin_0))
	endif
	if(iter_gw+iter_qp/=0) then
	  allocate(g_full(nbndf,nbndf,2,ndim3_tau,maxval(ndim3_k),
     &           nspin_0))
	  allocate(qp_spectr_weight(nbndf,maxval(ndim3_k),nspin_0))
      endif
	if(iter_gw+iter_qp/=0) then
	  allocate(sig_c_tau(nbndf,nbndf,2,ndim3_tau,maxval(ndim3_k),
     &           nspin_0))
	  allocate(sig_c_tau_old(nbndf,nbndf,2,ndim3_tau,
     &           maxval(ndim3_k),nspin_0))
        sig_c_tau_old=0.d0
      endif
	allocate(sig_c_0(nbndf,nbndf,maxval(ndim3_k),nspin_0))
	sig_c_0=0.d0
	allocate(z_ren(nbndf,nbndf,maxval(ndim3_k),nspin_0))
	if(iter_gw+iter_qp/=0) then
	  allocate(e_qp(nbndf,npnt,nspin_0))
	  allocate(z_qp(nfun,nbndf,maxval(ndim3_k),nspin_0))
	  allocate(q_qp(nbndf,nbndf,maxval(ndim3_k),nspin_0))
	endif
      if(iter_max_scf>iter_dft) then
	  allocate(v_h_bb_old(nbndf,nbndf,maxval(ndim3_k),nspin_0))
	  allocate(v_xc_bb(nbndf,nbndf,maxval(ndim3_k),nspin_0))
	  allocate(v_xc_bb_old(nbndf,nbndf,maxval(ndim3_k),nspin_0))
	  allocate(v_h_bb(nbndf,nbndf,maxval(ndim3_k),nspin_0))
	endif
	if(iter_hf+iter_gw+iter_qp/=0) then
	  allocate(sigx_solid_k(nbndf,nbndf,maxval(ndim3_k),nspin_0))
	  allocate(sigx_solid_k_old(nbndf,nbndf,maxval(ndim3_k),
     &           nspin_0))
        sigx_solid_k_old=0.d0
      endif
	if(magn.eq.2) then
	  allocate(spmt(maxmtb))
	  allocate(t_m(maxlfun,maxlfun,2,maxnsymb,nsort))
	  spmt=0.d0
      endif
	allocate(ntle_pb(0:maxpb,nsort))
      call timel('***** memory1_solid finished *******')
      end
