      subroutine sym_sig_r_from_nonsym_ii(f_in,f_out,ind_s)
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
      integer, intent(in) :: ind_s
	complex*16, intent(in) :: f_in(nr_full,max_r_star,nqdiv)
	complex*16, intent(out) :: f_out(nr_full,max_r_star,nqdiv)
	integer :: j,ig,i,ir,ir0,ind_j,ind_j0,i0,j0
      real*8 :: v(3),v1(3)
	f_out=(0.d0,0.d0)
	do ind_j=1,num_r_star(ind_s)
	  j=ind_r_star(ind_j,ind_s)
	  do i=1,nr_full
	    do ir=1,nqdiv
	      do ig=1,ngroup
	        v=rr_coarse(:,ir)+tshift_r(:,i,ig)-tshift_r(:,j,ig)
              call rotate(v(1),v(2),v(3),v1(1),v1(2),v1(3),u(2,ig),1)
	        call zone1_number(v1,gbas,ndiv,ir0)
	        i0=ip_r(i,ig)
	        j0=ip_r(j,ig)
	        ind_j0=ind_r_star_back(j0)
              f_out(i,ind_j,ir)=f_out(i,ind_j,ir)+f_in(i0,ind_j0,ir0)
	      enddo    !! over ig
	    enddo   !! over ir
	  enddo  !! over i
	enddo   !! over ind_j
	f_out=f_out/ngroup
      end