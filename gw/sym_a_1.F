      subroutine sym_a_1(ev_new,k,ev_bnd0,neig,ig,k0)
	use manager_mod
	use solid_mod
	use units_mod
      implicit none
      integer, intent(in) :: k,neig,ig,k0
	complex*16, intent(in) :: ev_bnd0(nbasmpw,neig)
	complex*16, intent(out) :: ev_new(nbasmpw,neig)
      integer :: j,ibnd,j0,k00,j1,gbs_number
      real*8 :: v(3),pi2,phase,gtild(3)
	complex*16 :: cc
	if(k==k0) then
	  ev_new=ev_bnd0
	  return
	endif
	pi2=pi+pi
      call rotate(pnt(1,k0),pnt(2,k0),pnt(3,k0),v(1),v(2),v(3),
     &            u(2,ig),2)
      gtild=pnt(:,k)-v    !! G_A
      k00=i_kref(k0)
	do j=1,nbask(k00)/nrel
	  j0=indgb(j,k)  !! G in GBS-list
	  v=gbs(:,j0)+gtild  !! G+G_A
        j1=gbs_number(v)   !! G+G_A in GBS-list
        j1=gbs_sym(j1,ig)   !! A^(-1)*(G+G_A) in GBS-list
        j1=iplf_bk(j1,k0)
	  v=pnt(:,k)+gbs(:,j0)  !! k+G
	  phase=-pi2*dot_product(v,shift(:,ig))
	  cc=dcmplx(cos(phase),sin(phase))
	  do ibnd=1,neig
	    ev_new(j,ibnd)=cc*ev_bnd0(j1,ibnd)
	  enddo
	enddo  !! over j
      end
