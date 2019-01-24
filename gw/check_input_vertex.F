      subroutine check_input_vertex
	use atom_mod
	use etot_mod
	use heg_mod
	use manager_mod
	use models_mod
	use parallel_mod
	use solid_mod
	use units_mod
	use vertex_mod
      implicit none
      if(iter_gwg==0) return
      if(mode_vertex_p/='0__ORDER'.and.w_sc_gwg=='non') then
        if(iter_qp/=0) then
          if(w_sc_qp=='scf') mode_vertex_p='BS____GW'
          if(w_sc_qp=='non') mode_vertex_p='BS___GW0'
        else if(iter_gw/=0) then
          if(w_sc_gw=='scf') mode_vertex_p='BS____GW'
          if(w_sc_gw=='non') mode_vertex_p='BS___GW0'
        else if(iter_hf/=0) then
          mode_vertex_p='BS____GV'
        else
          mode_vertex_p='BS___DFT'
        endif
      endif
      same_vertex=.true.
      if(mode_vertex_p/=mode_vertex_sigma) same_vertex=.false.
      if(iter_ladder_p/=iter_sigma_gwg) same_vertex=.false.
      full_vertex_p=.false.
      if(mode_vertex_p=='1__ORDER') full_vertex_p=.true.
      if(mode_vertex_p=='BS___GW0') full_vertex_p=.true.
      if(mode_vertex_p=='BS____GW') full_vertex_p=.true.
      if(mode_vertex_p=='BS__G3W2') full_vertex_p=.true.
      if(mode_vertex_p=='BS_G3W02') full_vertex_p=.true.
      full_vertex_sigma=.false.
      if(mode_vertex_sigma=='1__ORDER') full_vertex_sigma=.true.
      if(mode_vertex_sigma=='BS___GW0') full_vertex_sigma=.true.
      if(mode_vertex_sigma=='BS____GW') full_vertex_sigma=.true.
      if(mode_vertex_sigma=='BS__G3W2') full_vertex_sigma=.true.
      if(mode_vertex_sigma=='BS_G3W02') full_vertex_sigma=.true.
      if(mode_vertex_p=='1__ORDER') iter_ladder_p=1
      if(mode_vertex_sigma=='1__ORDER') iter_sigma_gwg=1
      end