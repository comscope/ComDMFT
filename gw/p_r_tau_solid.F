      subroutine p_r_tau_solid
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer :: ispin
	p_w_solid_q_tau=0.d0
	do ispin=1,nspin
	  call pol_0_spin_spin(ispin,ispin,p_w_solid_q_tau)
      enddo  !! over ispin
c ---------------------------------------------------------------      
	p_w_solid_q_tau=p_w_solid_q_tau*4.d0/nspin/nrel
      call timel('** P_R_TAU_SOLID finished **********')
      end