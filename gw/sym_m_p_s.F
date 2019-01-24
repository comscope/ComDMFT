	subroutine sym_m_p_s(xnew,xold,iq,iq0,ig)
	use atom_mod
	use manager_mod
	use solid_mod
	use units_mod
	implicit none
	integer, intent(in) :: iq,iq0,ig
	complex*16, intent(in) :: xold(n_pbtot)
	complex*16, intent(out) :: xnew(n_pbtot)
	integer :: i,i0,i1,gbs_number,iq00
	real*8 :: pi2,v(3),gtild(3),phase
	complex*16 :: cc
	iq00=i_kref(iq)
	pi2=pi+pi
      call rotate(pnt(1,iq0),pnt(2,iq0),pnt(3,iq0),v(1),v(2),v(3),
     &            u(2,ig),2)
      gtild=pnt(:,iq)-v
      call sym_mp_mt(ig,iq,xnew,xold,1)
      do i=1,nplwgw(iq00)
	  i0=indpw_gw(i,iq)
	  v=gbs(:,i0)+gtild  !! G+G_A
        i1=gbs_number(v)   !! G+G_A in GBS-list
        i1=gbs_sym(i1,ig)   !! A^(-1)*(G+G_A) in GBS-list
        i1=iplf_gk(i1,iq)
	  v=pnt(:,iq)+gbs(:,i0)
	  phase=-pi2*dot_product(v,shift(:,ig))
	  cc=dcmplx(cos(phase),sin(phase))
	  xnew(n_pbmt+i)=cc*xold(n_pbmt+i1)
	enddo
	end