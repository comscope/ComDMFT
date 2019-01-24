      subroutine sigma_x_mt_mt(nbas,ndim,sig_x0,yc,ind_y,pmm)
	use atom_mod
	use manager_mod
	use solid_mod
      implicit none
      complex*16, intent(in) :: pmm(nfun,nfun)
      integer :: ind_y(2*limlb,natom),nbas,ndim,iatom,isort,ndimb,
     &	       nt1,ind,ndimb2,indy,lm,l,lget,i,mj,li,ie,i0
      complex*16 :: sig_x0(ndim,ndim),yc(max_bs,nbas)
	complex*16, allocatable :: tmp(:,:)
	allocate(tmp(ndim,nfun))
	tmp=(0.d0,0.d0)
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
	      tmp(i,:)=pmm(i0,:)
	    enddo
	  enddo
        nt1=ind+nlctot(isort)+1
        indy=ind_y(1,iatom)
        ndimb2=2*ndimb
        call zgemm('c','n',nbas,nfun,ndimb2,(1.d0,0.d0),yc(indy,1),
     &	           max_bs,pmm(nt1,1),nfun,(1.d0,0.d0),
     &             tmp(nlcloc+1,1),ndim)
      enddo   !!! over iatom
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
	      sig_x0(:,i)=sig_x0(:,i)+tmp(:,i0)
	    enddo
	  enddo
        nt1=ind+nlctot(isort)+1
        indy=ind_y(1,iatom)
        ndimb2=2*ndimb
        call zgemm('n','n',ndim,nbas,ndimb2,(1.d0,0.d0),tmp(1,nt1),
     &	           ndim,yc(indy,1),max_bs,(1.d0,0.d0),
     &             sig_x0(1,nlcloc+1),ndim)
      enddo   !!! over iatom
	deallocate(tmp)
      end