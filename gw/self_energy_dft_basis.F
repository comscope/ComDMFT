      subroutine self_energy_dft_basis
	use atom_mod
	use etot_mod
      use manager_mod
      use parallel_mod
      use solid_mod
	use units_mod
	implicit none
	integer :: ispin,ind_k,k
c ------ Getting the matrix of Hartree potential  ----
  	  pv=0.d0
	  if(irel==2) pvj=0.d0
	  call vslli(0,1,v_mt_h,v_mt_xc,v_intr_h,v_intr_xc) ! only Hartree
	  do ispin=1,nspin
	    do ind_k=1,ndim3_k(me3_k+1)
	      k=n3_mpi_k(me3_k+1)+ind_k
	      call v_bloch(v_h_bb(1,1,ind_k,ispin),pnt(1,k),ispin,
     &                   n_bnd(k,ispin),nbndf,nbask(k),indgb(1,k),
     &                   ev_bnd(1,1,ind_k,ispin),z_bnd(1,1,ind_k,ispin),
     &                   pv,pvj,war,wari)
	    enddo
	  enddo
c ---------------------------------------------------------------	
        call output_static_band(v_h_bb,0)
c ------ Getting the exchange part of Self Energy  ----
        call sigx_k_solid
        call output_static_band(sigx_solid_k,2)
        call timel('**** SIGX_K_SOLID finished *********')
c ------ Getting the correlated part of Self Energy ----
	if(iter>iter_dft+iter_hf) then
          call sigc_r_tau_solid
          call timel('**** SIGC_R_TAU_SOLID finished *****')
      endif
c -------------------------------------------------------------------
      call timel('* SELF_ENERGY_DFT_BASIS finished ***')
      end
