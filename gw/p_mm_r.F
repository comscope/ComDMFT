      subroutine p_mm_r(p_mm,ind_tau,isp1,isp2)
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
      implicit none
      integer, intent(in) :: isp1,isp2,ind_tau
	real*8, intent(out) :: p_mm(n_pbmtm,n_pbmtm,nrr_max)
	integer :: isort,lm3,lm,l,ii,iii,indx,lget,icg1,icg2,icg,km3,lf3,
     &	         jatom,jsort,km2,lm2,lf2,l1m1,l1,km1,lm1,lf1,iatom,
     &           km4,lm4,lf4,ie,lf,km,ir0,ind_ir0
	real*8 :: cgg,eps,g1,g2,dum
	complex*16 :: g1c,g2c,cggc,dumc
	real*8, allocatable :: x(:,:),y(:),y1(:),gr1(:,:,:),gr2(:,:,:),
     &	                     p_mm_0(:,:)
	complex*16, allocatable :: xc(:,:),yc(:),gc1(:,:,:),gc2(:,:,:)
      if(irel<=1) then
        allocate(gr1(maxel,maxel,ndim_nrr(me3_k+1)))
        allocate(gr2(maxel,maxel,ndim_nrr(me3_k+1)))
	  call g_rs_from_ks_mm(gr1,dumc,1,ind_tau,isp1)
	  call g_rs_from_ks_mm(gr2,dumc,2,ind_tau,isp2)
      else if(irel==2) then
        allocate(gc1(maxel,maxel,ndim_nrr(me3_k+1)))
        allocate(gc2(maxel,maxel,ndim_nrr(me3_k+1)))
	  call g_rs_from_ks_mm(dum,gc1,1,ind_tau,isp1)
	  call g_rs_from_ks_mm(dum,gc2,2,ind_tau,isp2)
      endif	
	allocate(y1(n_pbmtm))
	if(irel<=1) then
	  allocate(x(maxel,n_pbmtm),y(n_pbmtm))
	  allocate(xc(1,1),yc(1))
	else if(irel==2) then
	  allocate(xc(maxel,n_pbmtm),yc(n_pbmtm))
	  allocate(x(1,1),y(1))
	endif
	eps=1.d-15
      p_mm=0.d0
      allocate(p_mm_0(n_pbmtm,n_pbmtm))
	do ind_ir0=1,ndim_nrr(me3_k+1)
	  ir0=n_mpi_nrr(me3_k+1)+ind_ir0
	  iatom=nrr_red(2,ir0)
	  jatom=nrr_red(3,ir0)
	  jsort=is(jatom)
	  isort=is(iatom)
	  p_mm_0=0.d0
	  do km3=1,lfunm(isort)
	    lm3=lm_isz(km3,isort)
	    lf3=lf_isz(km3,isort)
	    if(irel<=1) x=0.d0
	    if(irel==2) xc=(0.d0,0.d0)
	    do km4=1,lfunm(jsort)
            g2=0.0d0
            g2c=(0.0d0,0.0d0)
	      if(irel<=1) g2=gr2(km3,km4,ind_ir0)
	      if(irel==2) then
	        g2c=conjg(gc2(km3,km4,ind_ir0))
	        g2=abs(g2c)
	      endif
	      if(abs(g2)>=eps) then
	        lm4=lm_isz(km4,jsort)
	        lf4=lf_isz(km4,jsort)
	        do km2=1,lfunm(jsort)
	          lm2=lm_isz(km2,jsort)
	          lf2=lf_isz(km2,jsort)
                icg1=0
                icg2=0
	          if(irel<=1) then
	            ii = max0(lm4,lm2)
	            iii = min0(lm4,lm2)
                  indx = (ii*(ii-1))/2 + iii
                  icg1 = indxcg(indx)
                  icg2 = indxcg(indx+1) - 1
	          else
		          icg1=indxcgr(lm2,lm4,1,0)
                  icg2=indxcgr(lm2,lm4,2,0)
                endif 
                do icg = icg1, icg2
	            if(irel<=1) l1m1 = jcg(icg)
	            if(irel==2) l1m1 = jcgr(icg)
	            l1=lget(l1m1)
	            if(l1>lmpb(jsort)) cycle
	            if(irel<=1) then
	              cgg=g2*cg(icg)
	              do ie=1,ntle_pb(l1,jsort)
	                lf1=ind_prod0(ie,l1,jsort)
	                km1=ind_prod(ie,l1m1,jsort)
	                x(km2,km1)=x(km2,km1)+cgg
     &                          *fi0(lf2,lf4,lf1,jsort,isp1,isp2)
	              enddo
	            else if(irel==2) then
	              cggc=g2c*cgr(icg)
	              do ie=1,ntle_pb(l1,jsort)
	                lf1=ind_prod0(ie,l1,jsort)
	                km1=ind_prod(ie,l1m1,jsort)
	                xc(km2,km1)=xc(km2,km1)+cggc
     &                           *fi0(lf2,lf4,lf1,jsort,isp1,isp2)
	              enddo
	            endif
	          enddo
	        enddo  !! over km2
	      endif
          enddo   !! over km4
          do km1=1,lfunm(isort) 
	      lm1=lm_isz(km1,isort)
	      lf1=lf_isz(km1,isort)
              icg1=0
              icg2=0
	      if(irel<=1) then
	        y(1:n_pbmt0(jsort))=0.d0
	        do km2=1,lfunm(jsort)
	          g1=gr1(km1,km2,ind_ir0)
	          if(abs(g1)>=eps) then
	            y(1:n_pbmt0(jsort))=y(1:n_pbmt0(jsort))
     &		 	                     +g1*x(km2,1:n_pbmt0(jsort))
	          endif
              enddo   !! over km2
              ii = max0(lm3,lm1)
              iii = min0(lm3,lm1)
              indx = (ii*(ii-1))/2 + iii
              icg1 = indxcg(indx)
              icg2 = indxcg(indx+1) - 1
	      else if(irel==2) then
	        yc(1:n_pbmt0(jsort))=(0.d0,0.d0)
	        do km2=1,lfunm(jsort)
	          g1c=gc1(km1,km2,ind_ir0)
	          if(abs(g1c)>=eps) then
	            yc(1:n_pbmt0(jsort))=yc(1:n_pbmt0(jsort))
     &	  	 	                      +g1c*xc(km2,1:n_pbmt0(jsort))
	          endif
	        enddo   !! over km2
		      icg1=indxcgr(lm3,lm1,1,0)
              icg2=indxcgr(lm3,lm1,2,0)
            endif 
            do icg = icg1, icg2
	        if(irel<=1) lm = jcg(icg)
	        if(irel==2) lm = jcgr(icg)
	        l=lget(lm)
	        if(l>lmpb(isort)) cycle
	        do ie=1,ntle_pb(l,isort)
	          lf=ind_prod0(ie,l,isort)
	          km=ind_prod(ie,lm,isort)
	          if(irel<=1) then
	            cgg=fi0(lf1,lf3,lf,isort,isp1,isp2)*cg(icg)
		          y1(1:n_pbmt0(jsort))=y(1:n_pbmt0(jsort))*cgg
	          else if(irel==2) then
	            cggc=fi0(lf1,lf3,lf,isort,isp1,isp2)*cgr(icg)
		          y1(1:n_pbmt0(jsort))=yc(1:n_pbmt0(jsort))*cggc
		        endif
	          do ii=1,n_pbmt0(jsort)
	            p_mm_0(km,ii)=p_mm_0(km,ii)-y1(ii)
	          enddo
	        enddo
	      enddo  !! over icg
	    enddo  !! over km1
	  enddo  !! over km3
	  do iii=1,n_pbmt0(jsort)
	    do ii=1,n_pbmt0(isort)
	      p_mm(ii,iii,ir0)=p_mm_0(ii,iii)
	    enddo
	  enddo
	enddo  !! over ind_ir0
c ---------------------------------------------------------------
	deallocate(y1,p_mm_0,x,y,xc,yc)
	if(irel<=1) then
	  deallocate(gr1,gr2)
	else if(irel==2) then
	  deallocate(gc1,gc2)
	endif
	if(nproc_k/=1) call dgop(p_mm,n_pbmtm**2*nrr_max,'  +',comm_pnt)
      end
