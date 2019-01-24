      subroutine b_me
	use atom_mod
	use manager_mod
	use solid_mod
	use units_mod
      implicit none
	integer :: isort,iatom,nlmf,ind,j0,jfn,i0,ifn,iv,lm2,isymb,mtb,
     &	       irad,j,lmj,i,lmi,icg1,icg2,icg,ilm1_,ilm2_
      real*8 :: c2,cgg,cggj,dqdall
	real*8, allocatable :: t(:,:,:,:),work(:)
	allocate(t(maxlfun*(maxlfun+1)/2,maxnsymb,3,2),work(0:maxnrad))
	c2=clight**2
	do isort=1,nsort
        do iatom=1,natom
          if(is(iatom).eq.isort) goto 1
	  enddo
1	  nlmf=(lmpb(isort)+1)**2
        ind=0
	  do j0=1,lfun(isort)
	    jfn=ind_wf(j0,isort)
	    do i0=1,j0
	      ifn=ind_wf(i0,isort)
	      ind=ind+1
	      do iv=1,3
	        do lm2=1,nlmf
	          if(symb(lm2,iv,isort)) then
	            isymb=lmsymb(lm2,iv,isort)
                  mtb=indmtb(isymb,isort)
                  do irad=0,nrad(isort)
                    work(irad)=gfun(irad+ifn,1)*gfun(irad+jfn,1)
     &				          *b_mt(mtb+irad)*r(irad,isort)**2
     &                                       *dr(irad,isort)
	            enddo
                  t(ind,isymb,iv,1)=dqdall(h(isort),work,nrad(isort))
                  do irad=0,nrad(isort)
                    work(irad)=gfund(irad+ifn,1)*gfund(irad+jfn,1)
     &				        *b_mt(mtb+irad)*r(irad,isort)**2
     &                                     *dr(irad,isort)
	            enddo
                  t(ind,isymb,iv,2)=-dqdall(h(isort),work,nrad(isort))
     &				             /c2
	          endif
	        enddo   !! over lm2
	      enddo  !! over iv
	    enddo   !! over i0
	  enddo   !! over j0
	  do j=1,lfunm(isort)
	    j0=lf_isz(j,isort)
	    lmj=lm_isz(j,isort)
	    do i=1,j
	      i0=lf_isz(i,isort)
	      lmi=lm_isz(i,isort)
	      ind=j0*(j0-1)/2+i0
	      do iv=1,3
		    icg1=indxcgr(lmi,lmj,1,iv)
              icg2=indxcgr(lmi,lmj,2,iv)
              do icg=icg1,icg2
	          lm2=jcgr(icg)
	          cgg=dreal(cgr(icg))
	          cggj=dimag(cgr(icg))
                if(lm2.le.nlmf) then
	            if(symb(lm2,iv,isort)) then
	              isymb=lmsymb(lm2,iv,isort)
                    pv(i,j,iatom,1)=pv(i,j,iatom,1)
     &				               -cgg*t(ind,isymb,iv,1)
                    pvj(i,j,iatom)=pvj(i,j,iatom)-cggj*t(ind,isymb,iv,1)
	            endif
	          endif
	        enddo
	        ilm1_=lm_tild(lmi)
	        ilm2_=lm_tild(lmj)
		      icg1=indxcgr(ilm1_,ilm2_,1,iv)
              icg2=indxcgr(ilm1_,ilm2_,2,iv)
              do icg=icg1,icg2
	          lm2=jcgr(icg)
	          cgg=dreal(cgr(icg))
	          cggj=dimag(cgr(icg))
                if(lm2.le.nlmf) then
	            if(symb(lm2,iv,isort)) then
	              isymb=lmsymb(lm2,iv,isort)
                    pv(i,j,iatom,1)=pv(i,j,iatom,1)
     &				               -cgg*t(ind,isymb,iv,2)
                    pvj(i,j,iatom)=pvj(i,j,iatom)-cggj*t(ind,isymb,iv,2)
	            endif
	          endif
	        enddo  !! over icg
	      enddo   !! over iv
          enddo   !! over i
	  enddo   !! over j
	enddo  !!  over isort
	deallocate(t,work)
      end