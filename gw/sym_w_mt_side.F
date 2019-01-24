	subroutine sym_w_mt_side(key,ig,tr,w_0,n)
c     See identity VII-32-3
c     key = 1 - left transformation
c     key = 2 - right transformation
	use atom_mod
	use solid_mod
	use units_mod
	implicit none
	integer, intent(in) :: ig,key,n
	complex*16, intent(in) :: tr(natom)
	complex*16, intent(inout) :: w_0(n_pbmt,n)
	integer :: i,iatom,isort,i0
	real*8, allocatable :: sr(:,:),sj(:,:),s0(:,:)
	allocate(s0(n_pbmt,n),sj(n_pbmt,n),sr(n_pbmt,n))
	s0=real(w_0)
	call sym_pro_w_side(s0,sr,ip(1,ig),u(1,ig),n)
	s0=imag(w_0)
	call sym_pro_w_side(s0,sj,ip(1,ig),u(1,ig),n)
	w_0=dcmplx(sr,sj)
	deallocate(s0,sr,sj)
	do iatom=1,natom
	  isort=is(iatom)
	  i0=iopb(iatom)
	  i=i0+n_pbmt0(isort)-1
	  if(key==1) w_0(i0:i,:)=tr(iatom)*w_0(i0:i,:)
	  if(key==2) w_0(i0:i,:)=conjg(tr(iatom))*w_0(i0:i,:)
	enddo
	end