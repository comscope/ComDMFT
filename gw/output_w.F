      subroutine output_w(w_last)
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	logical, intent(in) :: w_last
	if(.not.maswrk) return
c ------- G=G'=0 components for q-points along K-line -----------
      call output_boson_q_nu(p_output,'__P')
      call output_boson_q_nu(w_output,'__W')
	if(w_last) then
        call output_boson_q_nu(e_output,'__E')
        call output_boson_q_nu(e1_output,'_E1')
      endif
111   format(a5,i4,1x,3f8.5,4x,a1)      
      end