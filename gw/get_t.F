      subroutine get_t
	use atom_mod
	use manager_mod
	use solid_mod
	use units_mod
      implicit none
	integer :: isort,iatom,nlmf,km,lm,lf,km1,lm1,lf1,ii,iii,indx,icg1,
     &	         icg2,icg,lm2,isym,isp,nlmb,iv,isymb,lm_,lm1_
      real*8 :: c2,cgg
      complex*16, allocatable :: gloc(:,:,:,:)
	allocate(gloc(maxel,maxel,natom,nspin_0))
	gloc=g_loc_0
	c2=clight**2
	t_ro=0.d0
	do isort=1,nsort
	  iatom=iat_1(isort)
	  nlmf=(lmpb(isort)+1)**2
        do km=1,lfunm(isort)
	    lm=lm_isz(km,isort)
	    lf=lf_isz(km,isort)
          do km1=1,lfunm(isort)
	      lm1=lm_isz(km1,isort)
	      lf1=lf_isz(km1,isort)
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
		        lm2 = jcg(icg)
	        else
	          lm2=jcgr(icg)
	        endif
	        if(lm2.le.nlmf) then
	          if(sym(lm2,isort)) then
	            isym=lmsym(lm2,isort)
	            do isp=1,nspin
                    cgg=dcmplx(0.0d0,0.0d0)
	              if(irel.le.1) then
	                cgg=cg(icg)*dreal(gloc(km,km1,iatom,isp))
	              else if(irel.eq.2) then
	                cgg=dreal(cgr(icg)*gloc(km,km1,iatom,isp))
	              endif
			        t_ro(lf,lf1,isym,isort,isp)=
     &			        t_ro(lf,lf1,isym,isort,isp)+cgg
	            enddo  !! over isp
	          endif
	        endif
	      enddo
	    enddo
	  enddo
	enddo
	t_ro=t_ro*2.d0/dfloat(nspin)/nrel
	if(magn.eq.2) then
	  t_m=0.d0
	  do isort=1,nsort
	    iatom=iat_1(isort)
	    nlmf=(lmpb(isort)+1)**2
	    nlmb=(lmb(isort)+1)**2
	    do iv=1,3
            do km=1,lfunm(isort)
	        lm=lm_isz(km,isort)
	        lf=lf_isz(km,isort)
              do km1=1,lfunm(isort)
	          lm1=lm_isz(km1,isort)
	          lf1=lf_isz(km1,isort)
c ----------- For P*P ------------------------------------------------
		        icg1=indxcgr(lm,lm1,1,iv)
                icg2=indxcgr(lm,lm1,2,iv)
                do icg = icg1, icg2
	            lm2=jcgr(icg)
	            if(lm2.le.nlmf.and.symb(lm2,iv,isort)) then
	              isymb=lmsymb(lm2,iv,isort)
	              cgg=dreal(cgr(icg)*gloc(km,km1,iatom,1))
			        t_m(lf,lf1,1,isymb,isort)=
     &			      t_m(lf,lf1,1,isymb,isort)+cgg
	            endif
	          enddo
c ----------- For Q*Q ------------------------------------------------
                lm_=lm_tild(lm)
	          lm1_=lm_tild(lm1)
	          if(lm_.gt.nlmb.or.lm1_.gt.nlmb) goto 3
		        icg1=indxcgr(lm_,lm1_,1,iv)
                icg2=indxcgr(lm_,lm1_,2,iv)
                do icg = icg1, icg2
	            lm2=jcgr(icg)
	            if(lm2.le.nlmf.and.symb(lm2,iv,isort)) then
	              isymb=lmsymb(lm2,iv,isort)
	              cgg=dreal(cgr(icg)*gloc(km,km1,iatom,1))
			        t_m(lf,lf1,2,isymb,isort)=
     &			    	t_m(lf,lf1,2,isymb,isort)-cgg/c2
	            endif
	          enddo
3			  continue
	        enddo
	      enddo
	    enddo
	  enddo
	endif
	deallocate(gloc)
      end
