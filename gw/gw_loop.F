      subroutine gw_loop
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      logical :: w_lda,w_old,w_rec
      integer :: it0,ierr,idum,n,ispin
      call set_adm(adm_sig)
c     ------------------------------------------------------------------ 
      it0=iter_dft+iter_hf
c     -----------------------------------------------------------------
      if(iter>it0) then
        w_lda=.false.
        w_old=.false.
!     ------------------------------------------------------------------
        allocate(p_w_solid_q_tau(n_pbtot,n_pbtot,ndim3_tau,
     &    maxval(ndim3_k)),stat=ierr)
        if(ierr/=0) call error_alloc('                       GW_LOOP',
     &    '               P_W_SOLID_Q_TAU',
     &    4,n_pbtot,n_pbtot,ndim3_tau,
     &    maxval(ndim3_k),idum,idum,
     &    idum,8)
c     -------------- GG - part -------------------------------------------
        call p_r_tau_solid
        call screened_w(w_lda)
        sig_c_tau=0.d0
      endif
c     --- Here we get GW part of Sigma in TAU, mix it -----------
c     ----------- and we get it on OMEGA (if need) -------      
      call self_energy
      if(iter>it0) deallocate(p_w_solid_q_tau)
c     ---- We get full Sigma on OMEGA -------     
      call self_energy_1
c     -----------------------------------------------------------------
c$$$  call search_mu_0(chem_pot,nelec,chem_pot)
      call search_mu_brent(chem_pot,nelec,chem_pot)            
      if(maswrk) write(iun,*)' Chemical Potential: ',chem_pot
      if(ubi==' qp') chem_pot_qp=chem_pot

      do ispin=1,nspin
        if(ubi==' hf') call output_eigen_states('__hf',2,ispin)
        if(ubi==' qp') call output_eigen_states('__qp',2,ispin)
        if(ubi==' gw') call output_eigen_states('__gw',2,ispin)
        if(ubi=='gwg') call output_eigen_states('_gwg',2,ispin)
      enddo      
      do ispin=1,nspin          
        if(ubi==' qp') call output_eigen_states('__qp',3,ispin)
      enddo          
      call timel('**** Search_mu finished ************')

! If doing dft or hf, then set_g calls green_0 and g_mt_from_g_full.  It also calls etot_gw_0, and etot_gw_1, and g_x_tau which is needed by etot_gw_1.  These calls to etot_* just calculate variables concerning the total energy.
      call set_g   

      call timel('******** Set_G finished ************')
      end
      
