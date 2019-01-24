      subroutine self_energy
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      call set_adm(adm_sig)
      if(ubi==' qp'.or.ubi==' gw') then	
        if(.not.allocated(sig_c_omega))
     &    allocate(sig_c_omega(nbndf,nbndf,2,ndim3_omega,
     &             maxval(ndim3_k),nspin))
        sig_c_omega=(0.d0,0.d0)
      endif
      call self_energy_dft_basis
c --------------- Mixing the GW part of Sigma --------------------
      call mixer_gw
c --------- Here we find a new spectrum ----------
      call bands_sig_0
c --------------------------------------------------------------
      v_h_bb_old=v_h_bb
      sigx_solid_k_old=sigx_solid_k
      if(ubi/=' hf') sig_c_tau_old=sig_c_tau
c -------------------------------------------------------------------
      call timel('**** SELF_ENERGY finished **********')
      end
