      subroutine gw_loop
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
	implicit none
	logical :: w_lda,w_old,w_rec
	integer :: it0,ierr,idum,n
      call set_adm(adm_sig)
c ------------------------------------------------------------------ 
      it0=iter_dft+iter_hf
c -----------------------------------------------------------------
      if(iter>it0) then
        w_lda=.false.
        w_old=.false.
! ------------------------------------------------------------------
	  allocate(p_w_solid_q_tau(n_pbtot,n_pbtot,ndim3_tau,
     &           maxval(ndim3_k)),stat=ierr)
        if(ierr/=0) call error_alloc('                       GW_LOOP',
     &                               '               P_W_SOLID_Q_TAU',
     &                               4,n_pbtot,n_pbtot,ndim3_tau,
     &                               maxval(ndim3_k),idum,idum,
     &                               idum,8)
c -------------- GG - part -------------------------------------------
        call p_r_tau_solid
        call screened_w(w_lda)
        sig_c_tau=0.d0
      endif
c --- Here we get GW part of Sigma in TAU, mix it -----------
c ----------- and we get it on OMEGA (if need) -------      
      call self_energy
      if(iter>it0) deallocate(p_w_solid_q_tau)
c ---- We get full Sigma on OMEGA -------     
      call self_energy_1
c -----------------------------------------------------------------
	call search_mu_0(chem_pot,nelec,chem_pot)
	if(maswrk) write(iun,*)' Chemical Potential: ',chem_pot
	if(ubi==' qp') chem_pot_qp=chem_pot
      call timel('**** Search_mu finished ************')
	call set_g    
      call timel('******** Set_G finished ************')
      end
