      subroutine integral_band_pw_band_vrt(isp1,isp2,ff,z_vrt,ev_vrt,
     &                                     b_pw_b,k,iq)
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer, intent(in) :: iq,k,isp1,isp2
	complex*16, intent(in) :: z_vrt(nfun,nbndf_bnd,2),
     &	                        ev_vrt(nbasmpw,nbndf_bnd,2),
     &                          ff(maxel,maxel,natom,nspin,nspin)
	complex*16, intent(out) :: b_pw_b(nbndf_bnd,nbndf_bnd)
      integer :: j,iatom,isort,k0,kq,kq0,ndimb,ind0,ind_kq,igb1,igb2,
     &           igb3,ig_tild,jbas,jnd,ia,ib,ic,ibas,ind,ia1,ib1,ic1
	real*8 :: v(3),g_tild(3),det,det3
	complex*16, allocatable :: s(:,:),tmp(:,:),tmp1(:,:)
	allocate(s(nbndtm,nbndtm))
	allocate(tmp(nbndtm,nbndf))
	allocate(tmp1(maxel,nbndf))
	b_pw_b=(0.d0,0.d0)
	k0=i_kref(k)
	v=pnt(:,k)-pnt(:,iq)
	call zone1_number(v,rb0,ndiv,ind_kq)
	kq=index_k1(ind_kq)
	kq0=i_kref(kq)
c ----------------- MT contribution ----------------------------------
	do iatom=1,natom
	  isort=is(iatom)
	  ndimb=lfunm(isort)
	  ind0=io_lem(iatom)
	  call zgemm('n','n',ndimb,n_low_bnd(k0,isp1),ndimb,
     &	           (1.d0,0.d0),ff(1,1,iatom,isp2,isp1),maxel,
     &             z_vrt(ind0,1,1),nfun,(0.d0,0.d0),tmp1,maxel)
	  call zgemm('c','n',n_low_bnd(kq0,isp2),
     &	           n_low_bnd(k0,isp1),ndimb,(1.d0,0.d0),
     &             z_vrt(ind0,1,2),nfun,tmp1,maxel,(1.d0,0.d0),
     &             b_pw_b,nbndf_bnd)
	enddo  !! over iatom
c ----------------- Interstitial contribution ------------------------
      det=det3(gbas(1,1),gbas(1,2),gbas(1,3))
	g_tild=pnt(:,k)-pnt(:,iq)-pnt(:,kq)
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
      call zgemm('n','n',nbask(kq0),n_low_bnd(k0,isp1),nbask(k0),
     &	             (1.d0,0.d0),s,nbndtm,ev_vrt(1,1,1),nbasmpw,
     &               (0.d0,0.d0),tmp,nbndtm)
	call zgemm('c','n',n_low_bnd(kq0,isp2),n_low_bnd(k0,isp1),
     &	         nbask(kq0),(1.d0,0.d0),ev_vrt(1,1,2),nbasmpw,tmp,
     &           nbndtm,(1.d0,0.d0),b_pw_b,nbndf_bnd)
	deallocate(s,tmp,tmp1)
      end