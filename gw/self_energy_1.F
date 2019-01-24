      subroutine self_energy_1
	use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
	use units_mod
	implicit none
	integer :: ispin,ind_k,k
	if(ubi==' qp'.or.ubi==' gw') then
        sig_c_omega=(0.d0,0.d0)
      endif
c --------- SIGMA_C (Omega) from SIGMA_C (Tau) -----------------------
	if(ubi==' qp'.or.ubi==' gw') then
        call sig_omega_solid
c -------------------------------------------------------------------
        call output_sigma_band(sig_c_omega,'sigc_w','sigc',evolt/2)
	endif
	if(ubi==' qp'.or.ubi==' gw') then
c ------- We prepare Sigma and its derivative at Omega=0 --------------
	  do ispin=1,nspin
	    do ind_k=1,ndim3_k(me3_k+1)
	      k=n3_mpi_k(me3_k+1)+ind_k
	      call sig_omega0_solid(n_bnd(k,ispin),ind_k,ispin)
	    enddo
	  enddo
      endif
c -------------------------------------------------------------------
      call timel('**** SELF_ENERGY_1 finished ********')
      end
