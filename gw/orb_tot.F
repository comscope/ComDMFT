      SUBROUTINE orb_tot(work,totmom,spmom,aorbmom,isort)
**********************************************************
* The calculation of the orbital and total moments       *
**********************************************************
	use atom_mod
	use manager_mod
	use solid_mod
	use units_mod
      IMPLICIT none
      INTEGER :: lf1,mt1,lf,mt,irad,isort,km1,lm1,l1,i1,mj1,li1,jj1,km,
     &	         lm,l,i,mj,li,jj,iatom
      real*8 :: totmom(3),work(0:maxnrad),spmom(3),aorbmom(3),ww,wwi,
     &	        f11,f22,c2,dqdall,spx,spy,spz,spxi,spyi,spzi,aorbx,
     &          aorby,aorbz,aorbxi,aorbyi,aorbzi,totx,toty,totz,totxi,
     &          totyi,totzi
	real*8, allocatable :: f1(:,:),f2(:,:)
      complex*16, allocatable :: gloc(:,:,:,:)
	allocate(gloc(maxel,maxel,natom,nspin_0))
	gloc=g_loc_0
	allocate(f1(maxlfun,maxlfun),f2(maxlfun,maxlfun))
	c2=clight**2
	spmom=0.d0
	aorbmom=0.d0
	totmom=0.d0
	iatom=iat_1(isort)
	do lf1=1,lfun(isort)
	  mt1=ind_wf(lf1,isort)
	  do lf=1,lfun(isort)
	    mt=ind_wf(lf,isort)
          do irad=0,nrad(isort)
            work(irad)=gfun(mt+irad,1)*gfun(mt1+irad,1)*dr(irad,isort)
     &		                                     *r(irad,isort)**2
          enddo
          f1(lf,lf1)=dqdall(h(isort),work,nrad(isort))
          do irad=0,nrad(isort)
            work(irad)=-gfund(mt+irad,1)*gfund(mt1+irad,1)
     &   	           *dr(irad,isort)*r(irad,isort)**2
          enddo
          f2(lf,lf1)=dqdall(h(isort),work,nrad(isort))/c2
        enddo  !!! over lf
      enddo  !!! over lf1
	do km1=1,lfunm(isort)
	  lm1=lm_isz(km1,isort)
	  lf1=lf_isz(km1,isort)
	  call getlimj(lm1,l1,i1,mj1,li1,0)
        jj1=l1+l1+i1
	  do km=1,lfunm(isort)
	    lm=lm_isz(km,isort)
	    lf=lf_isz(km,isort)
	    call getlimj(lm,l,i,mj,li,0)
          jj=l+l+i
          ww=dreal(gloc(km,km1,iatom,1))
          wwi=dimag(gloc(km,km1,iatom,1))
          f11=f1(lf,lf1)
          f22=f2(lf,lf1)
c -----------  The spin moment ---------------------------------------
          call sp_rel(l,i,mj,l1,i1,mj1,spx,spy,spz,spxi,spyi,spzi)
          spmom(1)=spmom(1)+f11*ww*spx
          spmom(2)=spmom(2)-f11*wwi*spyi
          spmom(3)=spmom(3)+f11*ww*spz
          call sp_rel(l+i,-i,mj,l1+i1,-i1,mj1,spx,spy,spz,spxi,spyi,
     &		        spzi)
          spmom(1)=spmom(1)+f22*ww*spx
          spmom(2)=spmom(2)-f22*wwi*spyi
          spmom(3)=spmom(3)+f22*ww*spz
c -----------  The orbital moment -------------------------------------
          call orb_rel(l,i,mj,l1,i1,mj1,aorbx,aorby,aorbz,aorbxi,aorbyi,
     &                 aorbzi)
          aorbmom(1)=aorbmom(1)+f11*ww*aorbx
          aorbmom(2)=aorbmom(2)-f11*wwi*aorbyi
          aorbmom(3)=aorbmom(3)+f11*ww*aorbz
          call orb_rel(l+i,-i,mj,l1+i1,-i1,mj1,aorbx,aorby,aorbz,aorbxi,
     &                 aorbyi,aorbzi)
          aorbmom(1)=aorbmom(1)+f22*ww*aorbx
          aorbmom(2)=aorbmom(2)-f22*wwi*aorbyi
          aorbmom(3)=aorbmom(3)+f22*ww*aorbz
c -----------  The total moment ---------------------------------------
	    call tot_rel(jj,mj,jj1,mj1,totx,toty,totz,totxi,totyi,totzi)
          totmom(1)=totmom(1)+f11*ww*totx
          totmom(2)=totmom(2)-f11*wwi*totyi
          totmom(3)=totmom(3)+f11*ww*totz
          totmom(1)=totmom(1)+f22*ww*totx
          totmom(2)=totmom(2)-f22*wwi*totyi
          totmom(3)=totmom(3)+f22*ww*totz
        enddo  !!! over km
      enddo  !!! over km1
	deallocate(f1,f2)
	deallocate(gloc)
      END