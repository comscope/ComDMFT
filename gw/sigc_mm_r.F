      subroutine sigc_mm_r(ispin,it,ind_tau,s_c)
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
      implicit none
	integer, intent(in) :: ispin,ind_tau,it
	complex*16, intent(out) :: s_c(nfun,nfun,ndim3_kk(me3_k+1))
	integer :: isort,lm3,lm,ii,iii,indx,lget,icg1,icg2,icg,km3,lf3,
     &	         jatom,jsort,njpb,km2,lm2,lf2,nipb,km1,lm1,lf1,ie,lf,km,
     &           ir0,iatom,lm22,l22,lf22,lm11,l11,lf11,km11,km22,isp,
     &           ind_ir0,npr,ndim0,ind_km1
	real*8 :: cgg,dum
	complex*16 :: cggc,dumc
	integer, allocatable :: ndim_km1(:),index_km1(:,:)
	real*8, allocatable :: x(:,:),y(:,:),sr(:,:,:),gr(:,:,:)
	complex*16, allocatable :: xc(:,:),yc(:,:),sc(:,:,:),wtmp(:,:),
     &	                         gc(:,:,:)
      isp=min(ispin,nspin_1)
c ---- Prepare Green's functions in real space ------------------------
      if(irel<=1) then
        allocate(gr(maxel,maxel,ndim_nrr(me3_k+1)))
	  call g_rs_from_ks_mm(gr,dumc,it,ind_tau,ispin)
      else if(irel==2) then
        allocate(gc(maxel,maxel,ndim_nrr(me3_k+1)))
	  call g_rs_from_ks_mm(dum,gc,it,ind_tau,ispin)
	endif
	if(irel<=1) then
	  allocate(sr(maxel,maxel,nrr_max))
	  sr=0.d0
	else if(irel==2) then
	  allocate(sc(maxel,maxel,nrr_max))
	  sc=(0.d0,0.d0)
	  allocate(wtmp(n_pbmtm,n_pbmtm))
	endif
	npr=1
	ndim0=maxel/npr
	if(ndim0*npr<maxel) ndim0=ndim0+1
	allocate(ndim_km1(npr))
	allocate(index_km1(ndim0,npr))
	do ind_ir0=1,ndim_nrr(me3_k+1)
	  ir0=n_mpi_nrr(me3_k+1)+ind_ir0
	  iatom=nrr_red(2,ir0)
	  jatom=nrr_red(3,ir0)
	  jsort=is(jatom)
	  isort=is(iatom)
	  nipb=n_pbmt0(isort)
	  njpb=n_pbmt0(jsort)
	  if(irel==2) wtmp=w_rs_mm(:,:,ind_ir0)
	  if(irel<=1) then
	    allocate(x(n_pbmtm,maxel),y(n_pbmtm,maxel))
	    allocate(xc(1,1),yc(1,1))
	  else if(irel==2) then
	    allocate(xc(n_pbmtm,maxel),yc(n_pbmtm,maxel))
	    allocate(x(1,1),y(1,1))
	  endif
	  call size_shift_par_load(lfunm(jsort),npr,ndim_km1,index_km1,
     &	                         ndim0)
        do ind_km1=1,ndim_km1(me_kk_pbr+1)
          km1=index_km1(ind_km1,me_kk_pbr+1)
	    lm1=lm_isz(km1,jsort)
	    lf1=lf_isz(km1,jsort)
	    if(irel<=1) x=0.d0
	    if(irel==2) xc=(0.d0,0.d0)
	    do km3=1,lfunm(jsort)
	      lm3=lm_isz(km3,jsort)
	      lf3=lf_isz(km3,jsort)
	      if(irel<=1) then
	        ii = max0(lm3,lm1)
	        iii = min0(lm3,lm1)
              indx = (ii*(ii-1))/2 + iii
              icg1 = indxcg(indx)
              icg2 = indxcg(indx+1) - 1 
	      else
		      icg1=indxcgr(lm3,lm1,1,0)
              icg2=indxcgr(lm3,lm1,2,0)
            endif
            do icg = icg1, icg2
	        if(irel<=1) lm22 = jcg(icg)
	        if(irel==2) lm22 = jcgr(icg)
	        l22=lget(lm22)
	        if(l22>lmpb(jsort)) cycle
	        do ie=1,ntle_pb(l22,jsort)
	          lf22=ind_prod0(ie,l22,jsort)
			    km22=ind_prod(ie,lm22,jsort)
	          if(irel<=1) then
	            cgg=fi0(lf3,lf1,lf22,jsort,isp,isp)*cg(icg)
	            x(km22,1:lfunm(isort))=x(km22,1:lfunm(isort))
     & 		        +gr(1:lfunm(isort),km3,ind_ir0)*cgg
	          else if(irel==2) then
	            cggc=fi0(lf3,lf1,lf22,jsort,isp,isp)*cgr(icg)
	            xc(km22,1:lfunm(isort))=xc(km22,1:lfunm(isort))
     & 		        +gc(1:lfunm(isort),km3,ind_ir0)*cggc
                endif
	        enddo   !! over ie
	      enddo   !! over icg
	    enddo   !! over km3
	    if(irel<=1) then
	      call dgemm('n','n',nipb,lfunm(isort),njpb,1.d0,
     &	 	           w_rs_mm(1,1,ind_ir0),n_pbmtm,x,n_pbmtm,0.d0,
     &                 y,n_pbmtm)
	    else if(irel==2) then
	      call zgemm('n','n',nipb,lfunm(isort),njpb,(1.d0,0.d0),wtmp,
     &	               n_pbmtm,xc,n_pbmtm,(0.d0,0.d0),yc,n_pbmtm)
          endif
	    do km=1,lfunm(isort)
	      lm=lm_isz(km,isort)
	      lf=lf_isz(km,isort)
	      do km2=1,lfunm(isort)
	        lm2=lm_isz(km2,isort)
	        lf2=lf_isz(km2,isort)
	        if(irel<=1) then
		        ii = max0(lm2,lm)
		        iii = min0(lm2,lm)
                indx = (ii*(ii-1))/2 + iii
                icg1 = indxcg(indx)
                icg2 = indxcg(indx+1) - 1
	        else
		        icg1=indxcgr(lm,lm2,1,0)
                icg2=indxcgr(lm,lm2,2,0)
              endif 
              do icg = icg1, icg2
	          if(irel<=1) lm11 = jcg(icg)
	          if(irel==2) lm11 = jcgr(icg)
	          l11=lget(lm11)
	          if(l11>lmpb(isort)) cycle
	          do ie=1,ntle_pb(l11,isort)
	            lf11=ind_prod0(ie,l11,isort)
	            km11=ind_prod(ie,lm11,isort)
	            if(irel<=1) then
	              cgg=fi0(lf2,lf,lf11,isort,isp,isp)*cg(icg)
	              sr(km,km1,ir0)=sr(km,km1,ir0)-cgg*y(km11,km2)
	            else if(irel==2) then
	              cggc=fi0(lf2,lf,lf11,isort,isp,isp)*cgr(icg)
	              sc(km,km1,ir0)=sc(km,km1,ir0)-cggc*yc(km11,km2)
	            endif
	          enddo
	        enddo
	      enddo
	    enddo  !! over km
	  enddo  !! over km1
        deallocate(x,y,xc,yc)
	enddo  !! over ind_ir0
	deallocate(ndim_km1,index_km1)
c --------------------------------------------------------------- 
	if(npr/=1) then
	  if(irel<=1) call dgop(sr,maxel**2*nrr_max,'  +',comm_tau_k)
	  if(irel==2) call dgop(sc,2*maxel**2*nrr_max,'  +',comm_tau_k)
	endif
	if(nproc_k/=1) then
	  if(irel<=1) call dgop(sr,maxel**2*nrr_max,'  +',comm_pnt)
	  if(irel==2) call dgop(sc,2*maxel**2*nrr_max,'  +',comm_pnt)
	endif
	if(irel<=1) then
	  deallocate(gr)
	else if(irel==2) then
	  deallocate(wtmp,gc)
	endif
	if(irel<=1) then
	  call sym_sig_rs(sr,s_c)
	  deallocate(sr)
	else if(irel==2) then
	  call sym_sig_rs_rel(sc,s_c)
	  deallocate(sc)
	endif
      end