      subroutine integral_band_pw_band(ff,znew,anew,znew1,anew1,b_pw_b,
     &                                 k,pn)
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer, intent(in) :: k
	real*8, intent(in) :: pn(3)
	complex*16, intent(in) :: znew(nfun,nbndf,nspin),
     &	                        anew(nbasmpw,nbndf,nspin),
     &                          znew1(nfun,nbndf,nspin),
     &                          anew1(nbasmpw,nbndf,nspin),
     &                          ff(maxel,maxel,natom,nspin,nspin)
	complex*16, intent(out) :: b_pw_b(nbndf,nbndf,nspin,nspin)
      integer :: j,iatom,isort,k0,kq,kq0,ndimb,ind0,ispin,jspin,
     &           igb1,igb2,igb3,ig_tild,jbas,jnd,ia,ib,ic,ibas,ind,ia1,
     &           ib1,ic1
	real*8 :: v(3),g_tild(3),det,det3
	complex*16, allocatable :: s(:,:),tmp(:,:),tmp1(:,:)
	allocate(s(nbndtm,nbndtm))
	allocate(tmp(nbndtm,nbndf))
	allocate(tmp1(maxel,nbndf))
	b_pw_b=(0.d0,0.d0)
	k0=i_kref(k)
	v=pnt(:,k)-pn
	call zone1_number(v,rb0,ndiv,kq)
	kq=index_k1(kq)
	kq0=i_kref(kq)
c ----------------- MT contribution ----------------------------------
      do jspin=1,nspin
        do ispin=1,nspin
	    do iatom=1,natom
	      isort=is(iatom)
	      ndimb=lfunm(isort)
	      ind0=io_lem(iatom)
	      call zgemm('n','n',ndimb,n_bnd(k0,jspin),ndimb,
     &	               (1.d0,0.d0),ff(1,1,iatom,ispin,jspin),maxel,
     &                 znew(ind0,1,jspin),nfun,(0.d0,0.d0),tmp1,maxel)
	      call zgemm('c','n',n_bnd(kq0,ispin),n_bnd(k0,jspin),ndimb,
     &	               (1.d0,0.d0),znew1(ind0,1,ispin),nfun,
     &                 tmp1,maxel,(1.d0,0.d0),b_pw_b(1,1,ispin,jspin),
     &                 nbndf)
	    enddo  !! over iatom
	  enddo  !! over ispin
	enddo   !! over jspin
c ----------------- Interstitial contribution ------------------------
      det=det3(gbas(1,1),gbas(1,2),gbas(1,3))
	g_tild=pnt(:,k)-pn-pnt(:,kq)
      igb1=idnint(det3(g_tild,gbas(1,2),gbas(1,3))/det)
      igb2=idnint(det3(gbas(1,1),g_tild,gbas(1,3))/det)
      igb3=idnint(det3(gbas(1,1),gbas(1,2),g_tild)/det)
	ig_tild=indplw(igb1,igb2,igb3)
c ------ We form the overlap matrix -----------------------------
	do jbas=1,nbask(k0)
        jnd=indgb(jbas,k)
        ia=igbs(1,ig_tild)+igbs(1,jnd)
        ib=igbs(2,ig_tild)+igbs(2,jnd)
        ic=igbs(3,ig_tild)+igbs(3,jnd)
        do ibas=1,nbask(kq0)
          ind=indgb(ibas,kq)
          ia1=ia-igbs(1,ind)
          ib1=ib-igbs(2,ind)
          ic1=ic-igbs(3,ind)
          j=indplw(ia1,ib1,ic1)
	    if(complex_ro) then
	      s(ibas,jbas)=dcmplx(sovr(j),sovi(j))
	    else
	      s(ibas,jbas)=dcmplx(sovr(j),0.d0)
	    endif
	  enddo
	enddo
	do jspin=1,nspin
	  do ispin=1,nspin
	    call zgemm('n','n',nbask(kq0),n_bnd(k0,jspin),nbask(k0),
     &	             (1.d0,0.d0),s,nbndtm,anew(1,1,jspin),nbasmpw,
     &               (0.d0,0.d0),tmp,nbndtm)
	    call zgemm('c','n',n_bnd(kq0,ispin),n_bnd(k0,jspin),
     &	             nbask(kq0),(1.d0,0.d0),anew1(1,1,ispin),nbasmpw,
     &               tmp,nbndtm,(1.d0,0.d0),b_pw_b(1,1,ispin,jspin),
     &               nbndf)
        enddo
      enddo   !!! over jspin
	deallocate(s,tmp,tmp1)
      end