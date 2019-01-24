      subroutine sym_r_r1_q_boson(vii1,vii0,iq,ig,ind_s)
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
      integer, intent(in) :: ind_s,iq,ig
	complex*16, intent(in) :: vii0(nr_full,max_r_star)
	complex*16, intent(out) :: vii1(nr_full,max_r_star)
	integer :: j0,i,j,j1,nsta,i1,iv(3)
	vii1=(0.d0,0.d0)
      nsta=num_r_star(ind_s)
	do j=1,nsta
	  j0=ind_r_star(j,ind_s)
	  j1=ip_r(j0,ig)
	  j1=ind_r_star_back(j1)
	  do i=1,nr_full
	    i1=ip_r(i,ig)
	    iv=itshift_r(:,i,ig)-itshift_r(:,j0,ig)
	    vii1(i,j)=phase_k_rr(iq,iv(1),iv(2),iv(3))*vii0(i1,j1)
	  enddo   !! over i
	enddo  !! over j
      end