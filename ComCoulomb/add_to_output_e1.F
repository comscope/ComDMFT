	subroutine add_to_output_e1(ind_k,k,n_pbt,e1,ind_nu)
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
	implicit none
	integer, intent(in) :: ind_k,k,n_pbt,ind_nu
	complex*16, intent(in) :: e1(n_pbt,n_pbt)
	integer :: i_nu
	complex*16 :: wint
	i_nu=me3_tau*ndim3_nu+ind_nu-1
c       ----------------------------------------------------------
	call pw_from_pb(ind_k,n_pbt,n_pbt,e1,wint,2)
	e1_output_int(i_nu,1,k)=dble(wint)
	e1_output_int(i_nu,2,k)=dble(e1(n_pbt,n_pbt))
	e1_output_int(i_nu,3,k)=dble(e1(1,1))
	e_output_int(i_nu,2,k)=1.d0/dble(wint)
	end
