      subroutine t_t1_x(isort,ispin,key)
c     key = 0 : for core levels
c     key = 1 : for valence levels      
	use atom_mod
	use manager_mod
	use solid_mod
	use units_mod
      implicit none
      integer, intent(in) :: key,isort,ispin
	integer :: nlm_c,lm1,li1,lm,li,ii,iii,indx,
     &	         icg1,icg2,icg3,icg4,icg,lm22,l22,km1,lf1,l1,i1,mj1,
     &           km2,lf2,li2,l2,i2,mj2,ii1,iii1,indx1,icgg,lm20,ic,l,
     &           lget,i,mj,lm2,iatom,lmc,lic,lc,mjc,nlm_all
	real*8 :: anl,df22,fac,facj,cgg,pi8
      complex*16, allocatable :: gloc(:,:,:)
	allocate(gloc(maxel,maxel,natom))
	gloc=g_loc_0(:,:,:,ispin)
	pi8=8.d0*pi
	nlm_c=nrel_core*(lc_max(isort)+1)**2
      nlm_all=0
	if(key==0) nlm_all=nlm_c
	if(key==1) nlm_all=nrel*(lmb(isort)+1)**2
	if(ncor(isort)/=0) then
	  t_x=0.d0
	  do lmc=1,nlm_c
	    if(irel_core.le.1) then
	      lic=lget(lmc)+1
	    else if(irel_core.eq.2) then
	      call getlimj(lmc,lc,ic,mjc,lic,0)
	    endif
	    do lm=1,nlm_all
              anl=0.0d0
	      if(irel_core.le.1) then
	        li=lget(lm)+1
	        anl=dfloat(li+li-1)
	      else if(irel_core.eq.2) then
	        call getlimj(lm,l,i,mj,li,0)
	        anl=dfloat(((li+1)/2)*2)
	      endif
	      if(irel_core.le.1) then
	        ii = max0(lm,lmc)
	        iii = min0(lm,lmc)
              indx = (ii*(ii-1))/2 + iii
              icg1 = indxcg(indx)
              icg2 = indxcg(indx+1) - 1 
	      else
	        icg1=indxcgr(lm,lmc,1,0)
              icg2=indxcgr(lm,lmc,2,0)
	      endif
            do icg = icg1, icg2
              if(irel_core.le.1) then
	          lm22=jcg(icg)
	          cgg=cg(icg)**2
	        else
	          lm22=jcgr(icg)
	          cgg=dconjg(cgr(icg))*cgr(icg)
	        endif
	        l22=lget(lm22)
	        if(key==0.and.l22>2*lc_max(isort)) cycle
	        if(key==1.and.l22>2*lmb(isort)) cycle
	        df22=-dfloat(l22+l22+1)*anl/pi8
		      t_x(li,lic,l22)=t_x(li,lic,l22)+cgg/df22
	      enddo
	    enddo      !! over lm - core lm
	  enddo   !! over lmc
	endif
	t1_x=0.d0
	iatom=iat_1(isort)
	do km1=1,lfunm(isort) 
	  lm1=lm_isz(km1,isort)
	  lf1=lf_isz(km1,isort)
	  if(irel.le.1) then
	    li1=lget(lm1)+1
	  else if(irel.eq.2) then
	    call getlimj(lm1,l1,i1,mj1,li1,0)
	  endif
	  do lm=1,nlm_all
          anl=0.0d0
	    if(irel_core.le.1) then
	      li=lget(lm)+1
	      anl=dfloat(li+li-1)
	    else if(irel_core.eq.2) then
	      call getlimj(lm,l,i,mj,li,0)
	      anl=dfloat(((li+1)/2)*2)
	    endif
	    if(irel.le.1) then
	      ii = max0(lm,lm1)
	      iii = min0(lm,lm1)
            indx = (ii*(ii-1))/2 + iii
            icg1 = indxcg(indx)
            icg2 = indxcg(indx+1) - 1 
	    else
	      icg1=indxcgr(lm,lm1,1,0)
            icg2=indxcgr(lm,lm1,2,0)
	    endif
          do icg = icg1, icg2
            if(irel.le.1) then
	        lm22=jcg(icg)
	        cgg=cg(icg)**2
	      else
	        lm22=jcgr(icg)
	        cgg=dconjg(cgr(icg))*cgr(icg)
	      endif
	      l22=lget(lm22)
	      df22=-dfloat(l22+l22+1)*anl/pi8
	      do km2=1,lfunm(isort)
	        lm2=lm_isz(km2,isort)
	        lf2=lf_isz(km2,isort)
	        if(irel.le.1) then
	          li2=lget(lm2)+1
	        else if(irel.eq.2) then
	          call getlimj(lm2,l2,i2,mj2,li2,0)
	        endif
	        if(irel.le.1) then
	          ii1 = max0(lm2,lm)
		        iii1 = min0(lm2,lm)
                indx1 = (ii1*(ii1-1))/2 + iii1
                icg3 = indxcg(indx1)
                icg4 = indxcg(indx1+1) - 1 
	        else
		        icg3=indxcgr(lm2,lm,1,0)
                icg4=indxcgr(lm2,lm,2,0)
	        endif
              do icgg = icg3, icg4
                if(irel.le.1) then
		        lm20 = jcg(icgg)
	          else
	            lm20=jcgr(icgg)
	          endif
	          if(lm20.eq.lm22) then
	            if(irel.le.1) then
	              fac=cg(icg)*cg(icgg)/df22
	            else
	              fac=dreal(cgr(icg)*cgr(icgg))/df22
	              facj=dimag(cgr(icg)*cgr(icgg))/df22
	            endif
	            if(irel.le.1) then
	              t1_x(l22,lf2,lf1,li)=t1_x(l22,lf2,lf1,li)
     &			    +dreal(gloc(km2,km1,iatom))*fac
	            else if(irel.eq.2) then
	              t1_x(l22,lf2,lf1,li)=t1_x(l22,lf2,lf1,li)
     &				+dreal(gloc(km2,km1,iatom)*dcmplx(fac,facj))
	            endif
	          endif
	        enddo
	      enddo    !! over km2
	    enddo
	  enddo      !! over lm - core lm
	enddo   !! over km1
	deallocate(gloc)
      end
