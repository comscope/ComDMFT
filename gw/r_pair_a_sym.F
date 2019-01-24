      subroutine r_pair_a_sym(ig,ir0,ir,iatom,jatom)
c
c	Gives the vector R=A^{-1}(R_0-T_t(A)+T_t'(A))
c
	use manager_mod
	use parallel_mod
      use sitegen_mod
	use solid_mod
      implicit none
	integer, intent(in) :: ig,ir0,iatom,jatom
	integer, intent(out) :: ir
      real*8 :: v(3),v1(3),v2(3)
      v=rr_pnt(1,ir0)*rbas(:,1)+rr_pnt(2,ir0)*rbas(:,2)
     &  		               +rr_pnt(3,ir0)*rbas(:,3)
	v1=v+tshift(:,iatom,ig)-tshift(:,jatom,ig)
      call rotate(v1(1),v1(2),v1(3),v2(1),v2(2),v2(3),u(2,ig),1)
	call zone1_number(v2,gbas,ndiv,ir)
      end