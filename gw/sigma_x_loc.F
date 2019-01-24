      subroutine sigma_x_loc(nbas,ndim,sig_x0,yc,ind_y)
	use atom_mod
	use manager_mod
	use solid_mod
      implicit none
      integer :: ind_y(2*limlb,natom),nbas,ndim,ibas0,iatom,
     &	         isort,ibas,nlct,ilc,jbas,jlc,ndimb,nt1,ind,ndimb2,nt
      complex*16 :: sig_x0(ndim,ndim),yc(max_bs,nbas)
	complex*16, allocatable :: tmp(:,:)
	sig_x0=(0.d0,0.d0)
c     &&& LO-LO contributions &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      ibas0=0
      do iatom=1,natom
        isort=is(iatom)
        ibas=ibas0
	  nlct=nlctot(isort)
        do ilc=1,nlct
          ibas=ibas+1
          jbas=ibas0
          do jlc=1,nlct
            jbas=jbas+1
	      sig_x0(ibas,jbas)=pv_x(ilc,jlc,iatom)
          enddo   !!! over jlc
        enddo   !!! over ilc
        ibas0=ibas0+nlct
      enddo   !!! over iatom
c     &&& PW-PW contributions &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	allocate(tmp(max_bs,nbas))
      do iatom=1,natom
        isort=is(iatom)
        ndimb=nrel*(lmb(isort)+1)**2
        nt1=nlctot(isort)+1
        ind=ind_y(1,iatom)
        ndimb2=2*ndimb
	  call zhemm('l','l',ndimb2,nbas,(1.d0,0.d0),pv_x(nt1,nt1,iatom),
     &	           maxel,yc(ind,1),max_bs,(0.d0,0.d0),tmp,max_bs)
        call zgemm('c','n',nbas,nbas,ndimb2,(1.d0,0.d0),yc(ind,1),
     &	           max_bs,tmp,max_bs,(1.d0,0.d0),
     &             sig_x0(nlcloc+1,nlcloc+1),ndim)
      enddo   !!! over iatom
	deallocate(tmp)
c     &&& PW-LO contributions &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      jbas=1
      do iatom=1,natom
        isort=is(iatom)
        ndimb=nrel*(lmb(isort)+1)**2
        nt=nlctot(isort)
        if(nt.ne.0) then
          nt1=nt+1
          ndimb2=2*ndimb
          ind=ind_y(1,iatom)
	    call zgemm('n','n',nt,nbas,ndimb2,(1.d0,0.d0),
     &		         pv_x(1,nt1,iatom),maxel,yc(ind,1),max_bs,
     &               (1.d0,0.d0),sig_x0(jbas,nlcloc+1),ndim)
          jbas=jbas+nt
        endif
      enddo   !!! over iatom
	sig_x0(nlcloc+1:ndim,1:nlcloc)=
     &	conjg(transpose(sig_x0(1:nlcloc,nlcloc+1:ndim)))
      end