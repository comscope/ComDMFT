      subroutine get_pc(ispin)
	use atom_mod
	use manager_mod
	use solid_mod
	use units_mod
      implicit none
	integer :: iatom,isort,lmm_c,l2,lic,lp1c,ic,kon,mtc,licmax,ispin,
     &	         ind,le1,mt1,ir,le,mt,km1,lm1,li1,lget,l1,i1,mj1,lmc,lc,
     &           mjc,ii,iii,indx,icg1,icg2,icg,lm2,km,lm,li,l,i,mj,ii1,
     &           iii1,indx1,icg3,icg4,icgg,lm22,idma,idmi,isp
      real*8 :: pi8,c2,allf,fac,ss,dqdall
	complex*16 :: cgg
	real*8, allocatable :: t(:,:,:),f(:),g(:),work(:),tmp(:),pp(:),
     &	                     qq(:)
      isp=min(ispin,nspin_1)
	licmax=1
	do isort=1,nsort
	  licmax=max(licmax,nrel*lc_max(isort)+1)
	enddo
	allocate(t(maxlfun*(maxlfun+1)/2,licmax,0:2*maxb),f(0:maxnrad),
     &	       g(0:maxnrad),work(0:maxnrad),tmp(0:maxnrad),
     &         pp(0:maxnrad),qq(0:maxnrad))
	pi8=8.d0*pi
	c2=clight**2
	allf=0.5d0*dfloat(nspin*nrel)
	do iatom=1,natom
	  isort=is(iatom)
	  lmm_c=nrel*(lc_max(isort)+1)**2   !! max LM for core
	  t=0.d0
	  do l2=0,2*lmb(isort)
	    fac=-pi8/dfloat(l2+l2+1)
	    ic=0
	    do lic=1,nrel*lc_max(isort)+1
	      if(irel.le.1) then
		      lp1c=lic
	      else
	        lp1c=lic/2+1
	      endif
            do kon=lp1c,konfig(lp1c-1,isort)-1
              ic=ic+1
	        mtc=indcor(ic,ispin,isort)
	        ind=0
	        do le1=1,lfun(isort)
	          mt1=ind_wf(le1,isort)
	          do ir=0,nrad(isort)
	            pp(ir)=gfun(mt1+ir,isp)*r(ir,isort)
	            if(irel.ge.1) qq(ir)=gfund(mt1+ir,isp)*r(ir,isort)
	          enddo
	          call var_int(dr(0,isort),pcor_old(mtc),qcor_old(mtc),
     &			             pp(0),qq(0),f,g,tmp(0),r(0,isort),
     &                       rl(0,l2,isort),work,nrad(isort)+1,h(isort))
	          do le=1,le1
	            ind=ind+1
	            mt=ind_wf(le,isort)
                  do ir=0,nrad(isort)
                    pp(ir)=gfun(mt+ir,isp)*pcor_old(mtc+ir)
                    if(irel.ge.1) pp(ir)=pp(ir)+gfund(mt+ir,isp)
     &			                               *qcor_old(mtc+ir)/c2
				  pp(ir)=pp(ir)*tmp(ir)*dr(ir,isort)*r(ir,isort)
                  enddo
                  ss=dqdall(h(isort),pp,nrad(isort))
	            t(ind,lic,l2)=t(ind,lic,l2)+fac*ss
	          enddo  !! over le
	        enddo  !! over le1
	      enddo   !! over kon
	    enddo    !! over lic
	  enddo   !! over l2
        do km1=1,lfunm(isort)
	    lm1=lm_isz(km1,isort)
	    if(irel.le.1) then
	      li1=lget(lm1)+1
	      l1=li1-1
	    else
	      call getlimj(lm1,l1,i1,mj1,li1,0)
	    endif
          do lmc=1,lmm_c
	      if(irel.le.1) then
		      lic=lget(lmc)+1
	        lc=lic-1
	      else
		      call getlimj(lmc,lc,ic,mjc,lic,0)
	      endif
	      if(irel.le.1) then
		      ii = max0(lmc,lm1)
		      iii = min0(lmc,lm1)
              indx = (ii*(ii-1))/2 + iii
              icg1 = indxcg(indx)
              icg2 = indxcg(indx+1) - 1 
	      else
		      icg1=indxcgr(lmc,lm1,1,0)
              icg2=indxcgr(lmc,lm1,2,0)
	      endif
            do icg = icg1, icg2
              if(irel.le.1) then
		        lm2 = jcg(icg)
	        else
	          lm2=jcgr(icg)
	        endif
	        l2=lget(lm2)
              do km=1,lfunm(isort)
	          lm=lm_isz(km,isort)
	          if(irel.le.1) then
	            li=lget(lm)+1
	            l=li-1
	          else
	            call getlimj(lm,l,i,mj,li,0)
	          endif
	          if(irel.le.1) then
		          ii1 = max0(lm,lmc)
		          iii1 = min0(lm,lmc)
                  indx1 = (ii1*(ii1-1))/2 + iii1
                  icg3 = indxcg(indx1)
                  icg4 = indxcg(indx1+1) - 1 
	          else
		          icg3=indxcgr(lm,lmc,1,0)
                  icg4=indxcgr(lm,lmc,2,0)
	          endif
                do icgg = icg3, icg4
                  if(irel.le.1) then
		            lm22 = jcg(icgg)
	            else
	              lm22=jcgr(icgg)
	            endif
	            if(lm22.eq.lm2) then
                    cgg=dcmplx(0.0d0,0.0d0)
	              if(irel.le.1) then
	                cgg=dcmplx(cg(icg)*cg(icgg),0.d0)
	              else if(irel.eq.2) then
	                cgg=cgr(icg)*cgr(icgg)
	              endif
	              le1=lf_isz(km1,isort)
	              le=lf_isz(km,isort)
	              idma=max(le1,le)
	              idmi=min(le1,le)
	              ind=idma*(idma-1)/2+idmi
	              pv_x(km,km1,iatom)=pv_x(km,km1,iatom)
     &		 		                  +cgg*t(ind,lic,l2)
	            endif
	          enddo  !! over icgg
	        enddo   !! over km
	      enddo    !! over icg
	    enddo    !! over kmc
	  enddo  !!! over km1
	enddo   !!! over iatom
	deallocate(t,f,g,work,tmp,pp,qq)
      end
