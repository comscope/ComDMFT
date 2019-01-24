      subroutine sym_tt_int(tt,n0,k,tt0,nop,ig,k0)
      use atom_mod
	use manager_mod
	use solid_mod
	use units_mod
      implicit none
      integer, intent(in) :: k,nop,ig,k0,n0
	complex*16, intent(in) :: tt0(n_pbtot,nop)
	complex*16, intent(out) :: tt(n0,nop)
      integer :: j,ibnd,j0,k00,j1,gbs_number
      real*8 :: pi2,v(3),gtild(3),phase
	complex*16 :: cc
	if(k==k0) then
	  tt=tt0(n_pbmt+1:n_pbmt+n0,:)
	  return
	endif
	pi2=pi+pi
      call rotate(pnt(1,k0),pnt(2,k0),pnt(3,k0),v(1),v(2),v(3),
     &            u(2,ig),2)
      gtild=pnt(:,k)-v    !! G_A
      k00=i_kref(k0)
	do j=1,nplwgw(k00)
	  j0=indpw_gw(j,k)  !! G in GBS-list
	  v=gbs(:,j0)+gtild  !! G+G_A
        j1=gbs_number(v)   !! G+G_A in GBS-list
        j1=gbs_sym(j1,ig)   !! A^(-1)*(G+G_A) in GBS-list
        j1=iplf_gk(j1,k0)
	  v=pnt(:,k)+gbs(:,j0)  !! k+G
	  phase=-pi2*dot_product(v,shift(:,ig))
	  cc=dcmplx(cos(phase),sin(phase))
	  do ibnd=1,nop
	    tt(j,ibnd)=cc*tt0(n_pbmt+j1,ibnd)
	  enddo
	enddo  !! over j
      end