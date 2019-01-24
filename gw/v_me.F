      subroutine v_me(ispin,l_st,key,v)
c	key=0 - only Hartree part
c	key=1 - Only XC part
c	key=2 - Hartree + XC
c  makes a nlm1*nlm2 matrix of p pertubation integrals.
c     l_st=1 : Start LM2 for V_HART calculation.
c     l_st=2 : Start LM2 for FLAPW matrix elements.
	use atom_mod
	use manager_mod
	use solid_mod
	use units_mod
      implicit none
	integer, intent(in) :: ispin,l_st,key
	real*8, intent(in) :: v(maxmt/nspin_0)
	integer :: isort,iatom,nlmf,ind,j0,jfn,i0,ifn,lm2,isym,mt,
     &	         irad,j,lmj,i,lmi,ix,icg1,icg2,icg,ima,imi,isp
      real*8 :: c2,cgg,cggj,dqdall,sqpi4,v0
	real*8, allocatable :: t(:,:),work(:)
      isp=min(ispin,nspin_1)
	allocate(t(maxlfun*(maxlfun+1)/2,maxnsym),work(0:maxnrad))
	c2=clight**2
	sqpi4=sqrt(4.d0*pi)
	do isort=1,nsort
	  iatom=iat_1(isort)
	  nlmf=(lmpb(isort)+1)**2
        ind=0
	  do j0=1,lfun(isort)
	    jfn=ind_wf(j0,isort)
	    do i0=1,j0
	      ifn=ind_wf(i0,isort)
	      ind=ind+1
	      do lm2=l_st,nlmf
	        if(sym(lm2,isort)) then
	          isym=lmsym(lm2,isort)
                mt=indmt(isym,isort,1)
                do irad=0,nrad(isort)
                  work(irad)=gfun(irad+ifn,isp)*gfun(irad+jfn,isp)
	          enddo
	          if(irel.ge.1) then
                  do irad=0,nrad(isort)
                    work(irad)=work(irad)+gfund(irad+ifn,isp)
     &			                     *gfund(irad+jfn,isp)/c2
	            enddo
	          endif
	          work(0)=0.d0
                do irad=1,nrad(isort)
	            v0=v(mt+irad)
	            if(lm2==1.and.key/=1.and.key/=3) v0=v0
     &				              -2.d0*sqpi4*z(isort)/r(irad,isort)
                  work(irad)=work(irad)*v0*dr(irad,isort)
     &				                    *r(irad,isort)**2
                enddo
                t(ind,isym)=dqdall(h(isort),work,nrad(isort))
	        endif
	      enddo   !! over lm2
	    enddo   !! over i0
	  enddo   !! over j0
	  do j=1,lfunm(isort)
	    j0=lf_isz(j,isort)
	    lmj=lm_isz(j,isort)
	    do i=1,j
	      i0=lf_isz(i,isort)
	      lmi=lm_isz(i,isort)
	      ima=max(i0,j0)
	      imi=min(i0,j0)
	      ind=ima*(ima-1)/2+imi
	      if(irel.le.1) then
              ix=max0(lmi,lmj)
              ix=(ix*(ix-1))/2+min0(lmi,lmj)
	        icg1=indxcg(ix)
	        icg2=indxcg(ix+1)-1
	      else
	        icg1=indxcgr(lmi,lmj,1,0)
              icg2=indxcgr(lmi,lmj,2,0)
	      endif
            do icg=icg1,icg2
              if(irel.le.1) then
		        lm2 = jcg(icg)
	          cgg=cg(icg)
	        else
	          lm2=jcgr(icg)
	          cgg=dreal(cgr(icg))
	          cggj=dimag(cgr(icg))
	        endif
              if(lm2.ge.l_st.and.lm2.le.nlmf) then
	          if(sym(lm2,isort)) then
	            isym=lmsym(lm2,isort)
                  pv(i,j,iatom,ispin)=pv(i,j,iatom,ispin)
     &				                 +cgg*t(ind,isym)
                  if(irel.eq.2) pvj(i,j,iatom)=pvj(i,j,iatom)
     &				                          +cggj*t(ind,isym)
	          endif
              endif
            enddo  !! over icg
	    enddo   !! over i
	  enddo   !! over j
	  do j=1,lfunm(isort)
	    do i=j+1,lfunm(isort)
            pv(i,j,iatom,ispin)=pv(j,i,iatom,ispin)
            if(irel.eq.2) pvj(i,j,iatom)=-pvj(j,i,iatom)
          enddo
	  enddo
	enddo  !! over isort
	deallocate(t,work)
      end