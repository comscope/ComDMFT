      subroutine sym_ii_rr_boson(vmi,tmp,ind_s)
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
      integer, intent(in) :: ind_s
	real*8, intent(inout) :: vmi(nr_full,max_r_star,nqdiv),
     &	                       tmp(nr_full,max_r_star,nqdiv)
	integer :: j0,i,ig,j,j1,ir,nsta,ir1,i1
	real*8 :: v(3),v1(3)
	tmp=0.d0
      nsta=num_r_star(ind_s)
	do ir=1,nqdiv
	  do ig=1,ngroup
	    do j=1,nsta
	      j0=ind_r_star(j,ind_s)
	      j1=ip_r(j0,ig)
	      j1=ind_r_star_back(j1)
	      do i=1,nr_full
	        i1=ip_r(i,ig)
	        v=rr_coarse(:,ir)+tshift_r(:,i,ig)-tshift_r(:,j0,ig)
	        call rotate(v(1),v(2),v(3),v1(1),v1(2),v1(3),u(2,ig),1)
              call zone1_number(v1,gbas,ndiv,ir1)
	        tmp(i,j,ir)=tmp(i,j,ir)+vmi(i1,j1,ir1)
	      enddo   !! over i
	    enddo  !! over j
	  enddo   !! over ig
	enddo   !! over ir
	vmi=tmp/ngroup
      end