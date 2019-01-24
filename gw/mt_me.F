      subroutine mt_me(ispin,a,ndim,ndim0,zbn,p,pj)
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
      implicit none
	integer, intent(in) :: ispin,ndim,ndim0
	real*8, intent(in) :: p(maxel,maxel,natom,nspin),
     &	                    pj(maxel,maxel,natom)
	complex*16, intent(in) :: zbn(nfun,nbndf)
	real*8, intent(out) :: a(ndim0,ndim0)
	integer :: iatom,isort,nl,lm1,lm,ind,ma
	complex*16, allocatable :: tmp(:,:),tm(:,:),tmp1(:,:)
	ma=maxel
	allocate(tmp(ma,ma),tm(ma,ndim))
	allocate(tmp1(ndim,ndim))
	tmp1=(0.d0,0.d0)
	do iatom=1,natom
	  isort=is(iatom)
	  nl=lfunm(isort)
	  if(irel.le.1) then
	    do lm1=1,nl
	      do lm=1,lm1
	        tmp(lm,lm1)=dcmplx(p(lm,lm1,iatom,ispin),0.d0)
	      enddo
	    enddo
	  else if(irel.eq.2) then
	    do lm1=1,nl
	      do lm=1,lm1
	        tmp(lm,lm1)=dcmplx(p(lm,lm1,iatom,ispin),
     &			                 pj(lm,lm1,iatom))
	      enddo
	    enddo
	  endif
c	tmp=conjg(tmp)
	  ind=io_lem(iatom)
	  call zhemm('L','U',nl,ndim,(1.d0,0.d0),tmp,ma,zbn(ind,1),nfun,
     &	           (0.d0,0.d0),tm,ma)
	  call zgemm('C','N',ndim,ndim,nl,(1.d0,0.d0),zbn(ind,1),nfun,tm,
     &	           ma,(1.d0,0.d0),tmp1,ndim)
	enddo
	call pack_hermit(tmp1,a,ndim,ndim,ndim0,0.d0,1.d0)
	deallocate(tmp,tm,tmp1)
      end