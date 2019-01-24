      subroutine sigma_x_mt_int(nbas,ndim,sig_x0,yc,ind_y,pmi)
	use atom_mod
	use manager_mod
	use solid_mod
      implicit none
      complex*16, intent(in) :: pmi(nfun,nbasmpw)
      integer :: ind_y(2*limlb,natom),nbas,ndim,iatom,isort,ndimb,
     &	         nt1,ind,ndimb2,indy,lm,l,lget,i,mj,li,ie,i0,j0,j
      complex*16 :: sig_x0(ndim,ndim),yc(max_bs,nbas)
	complex*16, allocatable :: sig_x1(:,:)
	allocate(sig_x1(nbas,nbas))
	sig_x1=(0.d0,0.d0)
      do iatom=1,natom
        isort=is(iatom)
	  ind=io_lem(iatom)-1
	  ndimb=nrel*(lmb(isort)+1)**2
	  do lm=1,ndimb
	    if(irel<=1) l=lget(lm)
	    if(irel==2) call getlimj(lm,l,i,mj,li,0)
	    do ie=1,ntle(l,isort)  !! over local functions
			if(augm(ie,l,isort)/='LOC') cycle
	      i=indbas(ie,lm,iatom)   !! position in hamiltonian matrix
	      i0=ind+indbasa(1,ie,lm,isort)   !! position in PVV
	      sig_x0(i,nlcloc+1:ndim)=sig_x0(i,nlcloc+1:ndim)
     &		                       +pmi(i0,1:nbas)
	    enddo
	  enddo
        nt1=ind+nlctot(isort)+1
        indy=ind_y(1,iatom)
        ndimb2=2*ndimb
        call zgemm('c','n',nbas,nbas,ndimb2,(1.d0,0.d0),yc(indy,1),
     &	           max_bs,pmi(nt1,1),nfun,(1.d0,0.d0),sig_x1,nbas)
      enddo   !!! over iatom
	sig_x0(nlcloc+1:ndim,1:nlcloc)=
     &	conjg(transpose(sig_x0(1:nlcloc,nlcloc+1:ndim)))
     	do j0=1,nbas
	  j=nlcloc+j0
	  do i0=1,j0
	    i=nlcloc+i0
	    sig_x0(i,j)=sig_x0(i,j)+sig_x1(i0,j0)+conjg(sig_x1(j0,i0))
	    sig_x0(j,i)=conjg(sig_x0(i,j))
	  enddo
	enddo
	deallocate(sig_x1)
      end