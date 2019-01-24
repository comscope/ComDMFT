	subroutine sym_w_mt_red(ig,pnt0,w_0,npbt)
c     See identity VII-32-3	
	use atom_mod
	use solid_mod
	use units_mod
	use vertex_mod
	implicit none
	integer, intent(in) :: npbt,ig
	complex*16, intent(inout) :: w_0(npbt,npbt)
	integer :: i,j,jatom,jsort,iatom,isort,j0,i0
	real*8 :: tt(3),pi2,fas,pnt0(3)
	complex*16 :: cf
	real*8, allocatable :: sr(:,:),sj(:,:),s0(:,:)
	allocate(s0(n_pbmt_red,n_pbmt_red),sj(n_pbmt_red,n_pbmt_red),
     &	       sr(n_pbmt_red,n_pbmt_red))
	do j=1,n_pbmt_red
	  do i=1,n_pbmt_red
	    s0(i,j)=dreal(w_0(i,j))
	  enddo
	enddo
	call sym_pro_w_red(s0,sr,ip(1,ig),u(1,ig))
	do j=1,n_pbmt_red
	  do i=1,n_pbmt_red
	    s0(i,j)=dimag(w_0(i,j))
	  enddo
	enddo
	call sym_pro_w_red(s0,sj,ip(1,ig),u(1,ig))
	do j=1,n_pbmt_red
	  do i=1,n_pbmt_red
	    w_0(i,j)=dcmplx(sr(i,j),sj(i,j))
	  enddo
	enddo
	deallocate(s0,sr,sj)
	pi2=pi+pi
	do jatom=1,natom
	  jsort=is(jatom)
	  j0=iopb_red(jatom)
	  j=j0+n_pbmt0_red(jsort)-1
	  do iatom=1,natom
	    isort=is(iatom)
	    i0=iopb_red(iatom)
	    i=i0+n_pbmt0_red(isort)-1
	    tt=tshift(:,iatom,ig)-tshift(:,jatom,ig)
	    fas=pi2*dot_product(pnt0,tt)
	    cf=dcmplx(cos(fas),sin(fas))
	    w_0(i0:i,j0:j)=cf*w_0(i0:i,j0:j)
	  enddo
	enddo
	end