      subroutine get_pv(ispin)
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
      integer, intent(in) :: ispin
	integer :: iatom,isort,l2,ind1,le4,mt4,ir,le1,mt1,ind,le,mt,
     &	         le3,mt3,km4,lm4,lf4,km1,lm1,lf1,idma,idmi,ii,iii,indx,
     &           icg1,icg2,icg,lm2,lget,km,lm,lf,km3,lm3,lf3,ii1,iii1,
     &           indx1,icg11,icg22,icgg,lm22,l,ie,jn,ind_km4,ind_km1,
     &           isp
	real*8 :: pi8,c2,fac,ss,dqdall,err
	complex*16 :: cgg
      integer, allocatable :: nd_k(:),nm_k(:),nd_t(:),nm_t(:)
	real*8, allocatable :: t(:,:,:),f(:),g(:),work(:),tmp(:),pp(:),
     &	                     qq(:),pp1(:),qq1(:)
      complex*16, allocatable :: tm(:,:),gloc(:,:,:)
	allocate(gloc(maxel,maxel,natom))
	gloc=g_loc_0(:,:,:,ispin)
      isp=min(ispin,nspin_1)
	allocate(t(maxlfun*(maxlfun+1)/2,maxlfun*(maxlfun+1)/2,0:2*maxb),
     &	       f(0:maxnrad),g(0:maxnrad),work(0:maxnrad),tmp(0:maxnrad),
     &         pp(0:maxnrad),qq(0:maxnrad),pp1(0:maxnrad),
     &         qq1(0:maxnrad))
	pi8=8.d0*pi
	c2=clight**2
	allocate(nd_k(nproc_k))
	allocate(nm_k(nproc_k))
	allocate(nd_t(nproc_tau))
	allocate(nm_t(nproc_tau))
	do iatom=1,natom
	  isort=is(iatom)
	  t=0.d0
	  do l2=0,2*lmb(isort)
	    fac=-pi8/dfloat(l2+l2+1)
	    ind1=0
	    do le4=1,lfun(isort)
	      mt4=ind_wf(le4,isort)
	      do ir=0,nrad(isort)
	        pp(ir)=gfun(mt4+ir,isp)*r(ir,isort)
	        if(irel.ge.1) qq(ir)=gfund(mt4+ir,isp)*r(ir,isort)
	      enddo
	      do le1=1,le4
	        ind1=ind1+1
	        mt1=ind_wf(le1,isort)
              do ir=0,nrad(isort)
                pp1(ir)=gfun_old(mt1+ir,isp)*r(ir,isort)
                if(irel.ge.1)qq1(ir)=gfund_old(mt1+ir,isp)*r(ir,isort)
              enddo
	        call var_int(dr(0,isort),pp,qq,pp1,qq1,f,g,tmp,r(0,isort),
     &			           rl(0,l2,isort),work(0),nrad(isort)+1,
     &                     h(isort))
	        ind=0
	        do le=1,lfun(isort)
	          mt=ind_wf(le,isort)
	          do le3=1,le
	            ind=ind+1
	            mt3=ind_wf(le3,isort)
                  do ir=0,nrad(isort)
                    work(ir)=gfun(mt+ir,isp)*gfun_old(mt3+ir,isp)
                    if(irel.ge.1) work(ir)=work(ir)+gfund(mt+ir,isp)
     &                                    *gfund_old(mt3+ir,isp)/c2
  				  work(ir)=work(ir)*tmp(ir)*dr(ir,isort)
     &				                       *r(ir,isort)**2
                  enddo
                  ss=dqdall(h(isort),work,nrad(isort))
	            t(ind,ind1,l2)=fac*ss
	          enddo
	        enddo
	      enddo
	    enddo
	    call check_sym_matr(t(1,1,l2),lfun(1)*(lfun(1)+1)/2,
     &	    maxlfun*(maxlfun+1)/2,err)
	  enddo   !! over l2
	  call size_shift_par(lfunm(isort),nproc_k,nd_k,nm_k)
	  call size_shift_par(lfunm(isort),nproc_tau,nd_t,nm_t)
	  allocate(tm(lfunm(isort),lfunm(isort)))
	  tm=(0.d0,0.d0)
	  do ind_km4=1,nd_t(me3_tau+1)
	    km4=nm_t(me3_tau+1)+ind_km4
	    lm4=lm_isz(km4,isort)
	    lf4=lf_isz(km4,isort)
	    do ind_km1=1,nd_k(me3_k+1)
	      km1=nm_k(me3_k+1)+ind_km1
	      lm1=lm_isz(km1,isort)
	      lf1=lf_isz(km1,isort)
	      idma=max(lf4,lf1)
	      idmi=min(lf4,lf1)
	      ind1=idma*(idma-1)/2+idmi
	      if(irel.le.1) then
		      ii = max0(lm1,lm4)
		      iii = min0(lm1,lm4)
              indx = (ii*(ii-1))/2 + iii
              icg1 = indxcg(indx)
              icg2 = indxcg(indx+1) - 1 
	      else
		      icg1=indxcgr(lm1,lm4,1,0)
              icg2=indxcgr(lm1,lm4,2,0)
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
	          lf=lf_isz(km,isort)
                do km3=1,lfunm(isort)
	            lm3=lm_isz(km3,isort)
	            lf3=lf_isz(km3,isort)
	            idma=max(lf3,lf)
	            idmi=min(lf3,lf)
	            ind=idma*(idma-1)/2+idmi
	            if(irel.le.1) then
		            ii1 = max0(lm3,lm)
		            iii1 = min0(lm3,lm)
                    indx1 = (ii1*(ii1-1))/2 + iii1
                    icg11 = indxcg(indx1)
                    icg22 = indxcg(indx1+1) - 1 
	            else
		            icg11=indxcgr(lm3,lm,1,0)
                    icg22=indxcgr(lm3,lm,2,0)
	            endif
                  do icgg = icg11, icg22
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
	                cgg=cgg*gloc(km1,km,iatom)
	                tm(km3,km4)=tm(km3,km4)+cgg*t(ind,ind1,l2)
	              endif
	            enddo  !! over icgg
	          enddo   !! over km3
	        enddo   !! over km
	      enddo   !! over icg
	    enddo   !! over km1
	  enddo   !! over km4
	  if(nproc_k/=1) call dgop(tm,2*lfunm(isort)**2,'  +',comm_pnt)
	  if(nproc_tau/=1) call dgop(tm,2*lfunm(isort)**2,'  +',comm_pnk)
	  pv_x(1:lfunm(isort),1:lfunm(isort),iatom)=
     &	  pv_x(1:lfunm(isort),1:lfunm(isort),iatom)+tm
	  deallocate(tm)
	enddo   !!! over iatom
c ----- The initial shift for ferromagnetic Cerium ---------------
      if(irel<2.and.natom==1.and.txtel(1)=='Ce  ') then
        if(nspin==2.and.abs(b_extval)<1.d-8) then
          if(iter-iter_dft<=3) then
            lm=11
            l=lget(lm)
            do ie=1,ntle(l,1)
              if(augm(ie,l,1)=='APW') then
                do jn=1,2
	            km=indbasa(jn,ie,lm,1)
	            ss=magn_shift(1)*ffsmt(jn,jn,ie,ie,l+1,1,1)
	            ss=ss*(-1)**ispin
                  pv_x(km,km,1)=pv_x(km,km,1)-ss
                enddo
              endif
            enddo
          endif
        endif
      endif	
	deallocate(t,f,g,work,tmp,pp,qq,pp1,qq1)
	deallocate(nd_k,nm_k,nd_t,nm_t,gloc)
      end
