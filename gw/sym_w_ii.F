      subroutine sym_w_ii(iq,ig,w,n,n0)     
	use atom_mod
	use manager_mod
	use parallel_mod
	use units_mod
	use solid_mod
      implicit none
	integer, intent(in) :: n,n0,iq
	complex*16, intent(inout)  :: w(n0,n0)
	integer :: ig,j,i,j0,j1,i0,i1,gbs_number,igtild,iq0
	real*8 :: pi2,v(3),phase
	complex*16 :: cc
	complex*16, allocatable :: tmp(:,:)
	pi2=pi+pi
	allocate(tmp(n,n))
	iq0=k_sym_0(iq,ig)
	igtild=g_sym_0(iq,ig)
	tmp=w(1:n,1:n)
	w=(0.d0,0.d0)
      do j=1,n
	  j0=indpw_gw(j,iq)  !! G' in GBS-list
	  v=gbs(:,j0)+gbs(:,igtild)  !! G'+G_A
        j1=gbs_number(v)   !! G'+G_A in GBS-list
        j1=gbs_sym(j1,ig)   !! A^(-1)*(G'+G_A) in GBS-list
        j1=iplf_gk(j1,iq0)
        do i=1,n
	    i0=indpw_gw(i,iq)
	    v=gbs(:,i0)+gbs(:,igtild)  !! G+G_A
          i1=gbs_number(v)   !! G+G_A in GBS-list
          i1=gbs_sym(i1,ig)   !! A^(-1)*(G+G_A) in GBS-list
          i1=iplf_gk(i1,iq0)
	    v=gbs(:,j0)-gbs(:,i0)
	    phase=pi2*dot_product(v,shift(:,ig))
	    cc=dcmplx(cos(phase),sin(phase))
	    w(i,j)=cc*tmp(i1,j1)
	  enddo
	enddo   !! over j
	deallocate(tmp)
      end
