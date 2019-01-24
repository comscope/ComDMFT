	subroutine prod_b2(ndim0,key,nl,ind_l,limax)
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
c	key=0 - Run to find the dimensions
c	key=1 - Run to fill the arrays
      integer, intent(in) :: ndim0,limax,nl(limax,limax),
     &                       ind_l(limax,limax,maxpb+1),key
	integer :: ii,isort,n1,n2,l2,n0,ispin,in,li,mt,ir,l,jn,li1,mt1,l1,
     &	         nl0,j,i,lwork,il,iu,nev,info,n,i0,m2,lm2,ii2,iatom,
     &           jspin
	real*8 :: c2,s00,dqdall,abstol,vl,vu,sqe
	integer, allocatable :: iwork(:),ifail(:)
	real*8, allocatable :: s0(:,:),evec(:,:),eval(:),pp(:),pp1(:),
     &	                     qq(:),qq1(:),pb(:,:),work(:),wrk(:)
	allocate(s0(ndim0,ndim0),evec(ndim0,ndim0),eval(ndim0),
     &	       pp(0:maxnrad),pp1(0:maxnrad),qq(0:maxnrad),
     &         qq1(0:maxnrad),pb(0:maxnrad,ndim0),work(0:maxnrad))
	c2=clight**2
        n1=0
        n2=0
	if(key.eq.0) then
	  n_pbmt=0
	  maxpbr=0
	  max_e_pb=0
	endif 
	if(key==1) fun_pb=0.d0
	ii=1
	do isort=1,nsort
	  if(key==0.and.maswrk) write(iun,*)
     &	  ' Product basis for isort ',isort
	  if(key.eq.0) then
	    lfun_pb(isort)=0
	    n_pbmt0(isort)=0
	  endif   
	  if(key.eq.1) then
	    n1=0
	    n2=0
	  endif
	  do l2=0,lmpb(isort)
	    if(key==0.and.maswrk) write(iun,*)' Eigen values for L = ',l2
	    n0=0
	    do ispin=1,nspin
c	      do in=1,ncor(isort)+lfun(isort)
c	        if(in.le.ncor(isort)) then
c			    li=li_cor(in,isort)
c	          mt=indcor(in,ispin,isort)
c	          do ir=1,nrad(isort)
c	            pp(ir)=pcor(mt+ir)/r(ir,isort)
c	            if(irel.ge.1) qq(ir)=qcor(mt+ir)/r(ir,isort)
c	          enddo
c		        call fit_zero(pp,r(0,isort))
c		        if(irel>=1) call fit_zero(qq,r(0,isort))
c	        else if(in.gt.ncor(isort)) then
c			    li=li_val(in-ncor(isort),isort)
c	          mt=ind_wf(in-ncor(isort),isort)
c	          do ir=0,nrad(isort)
c	            pp(ir)=gfun(mt+ir,ispin)
c	            if(irel.ge.1) qq(ir)=gfund(mt+ir,ispin)
c	          enddo
c	        endif	
	        
	      do in=1,lfun(isort)
 			  li=li_val(in,isort)
	        mt=ind_wf(in,isort)
	        do ir=0,nrad(isort)
	          pp(ir)=gfun(mt+ir,ispin)
	          if(irel.ge.1) qq(ir)=gfund(mt+ir,ispin)
	        enddo

	        if(irel.le.1) l=li-1
	        if(irel.eq.2) l=li/2
	      do jspin=1,ispin
	        do jn=1,lfun(isort)
			    li1=li_val(jn,isort)
	          mt1=ind_wf(jn,isort)
	          do ir=0,nrad(isort)
	            pp1(ir)=gfun_old(mt1+ir,jspin)
	            if(irel.ge.1) qq1(ir)=gfund_old(mt1+ir,jspin)
	          enddo
	          if(irel.le.1) l1=li1-1
	          if(irel.eq.2) l1=li1/2
	          do nl0=1,nl(li,li1)
	            if(ind_l(li,li1,nl0).eq.l2) then
	              n0=n0+1
	              do ir=0,nrad(isort)
	                pb(ir,n0)=pp(ir)*pp1(ir)
	                if(irel.ge.1) pb(ir,n0)=pb(ir,n0)
     &					                     +qq(ir)*qq1(ir)/c2
	              enddo
                    do ir=0,nrad(isort)
                      work(ir)=pb(ir,n0)*pb(ir,n0)*dr(ir,isort)
     &						                      *r(ir,isort)**2
                    enddo
                    s00=dqdall(h(isort),work,nrad(isort))
	              s00=1.d0/sqrt(s00)
	              call dscal(nrad(isort)+1,s00,pb(0,n0),1)
	            endif
	          enddo
	        enddo   !! over jn
	    enddo   !! over jspin
	      enddo   !! over in
	    enddo   !! over ispin
	    do j=1,n0
	      do i=1,j
              do ir=0,nrad(isort)
                pp(ir)=pb(ir,i)*pb(ir,j)*dr(ir,isort)*r(ir,isort)**2
              enddo
              s0(i,j)=dqdall(h(isort),pp,nrad(isort))
	      enddo
	    enddo
c ---------- Diagonalize ----------------------------------------
	    abstol=1.d-70
          lwork=max(1,8*n0)
          allocate(wrk(lwork),iwork(5*n0),ifail(n0))
	    call DSYEVX('V','A','U',n0,s0,ndim0,VL,VU,IL,IU,ABSTOL,nev,
     &		          eval,evec,ndim0,wrk,lwork,iwork,ifail,INFO) 
	    if(info.ne.0) then
	      if(maswrk) write(iun,*)'Prod_b2: Info = ',info
	      if(info.gt.n0) call ending
	    endif
          deallocate(wrk,iwork,ifail)
	    n=0
	    do i=n0,1,-1
	      if(eval(i)/eval(n0).ge.eps_pb) n=n+1
	      if(n==lim_pb_mt(l2,isort)) exit
	    enddo
	    if(key.eq.0) then
	      ntle_pb(l2,isort)=n
	      lfun_pb(isort)=lfun_pb(isort)+n
	      n_pbmt0(isort)=n_pbmt0(isort)+n*(l2+l2+1)
	      max_e_pb=max(n,max_e_pb)
	      if(maswrk) then
	        write(iun,'(7(1x,f9.6))')(eval(i)/eval(n0),i=n0,n0-n+1,-1)
	      endif
	    else if(key.eq.1) then
	      do i=1,n  
	        n1=n1+1
	        ind_pb(n1,isort)=ii
	        lval_pb(n1,isort)=l2
	        ind_prod0(i,l2,isort)=n1
	        i0=n0-i+1
	        sqe=sqrt(eval(i0))
	        do j=1,n0
	          do ir=0,nrad(isort)
	            fun_pb(ii+ir)=fun_pb(ii+ir)+pb(ir,j)*evec(j,i0)/sqe
	          enddo
	        enddo
	        ii=ii+nrad(isort)+1
	        do m2=-l2,l2
	          n2=n2+1
	          lm2=l2*(l2+1)+m2+1
	          lm_pbmt(n2,isort)=lm2
	          lf_pb(n2,isort)=n1
	          ind_prod(i,lm2,isort)=n2
	        enddo
	      enddo
	    endif
	  enddo   !! over l2
	  if(key.eq.0) then
	    maxlfpb=max(maxlfpb,lfun_pb(isort))
	    n_pbmtm=max(n_pbmtm,n_pbmt0(isort))
	    n_pbmt=n_pbmt+n_pbmt0(isort)*nhsort(isort)
	    maxpbr=maxpbr+lfun_pb(isort)*(nrad(isort)+1)
	  endif
	enddo   !! over isort
	if(key.eq.0) then
	  ii2=1
	  do iatom=1,natom
	    isort=is(iatom)
	    iopb(iatom)=ii2
	    ii2=ii2+n_pbmt0(isort)
	  enddo
	  ii2=1
	  do isort=1,nsort
	    iopb0(isort)=ii2
	    ii2=ii2+n_pbmt0(isort)
	  enddo
	  if(maswrk) then
	    do isort=1,nsort
	      write(iun,*)' Product basis info for ',txtel(isort)
	      write(iun,'(1x,a1,4x,a4,1x,a11)')'L','ntle','ntle*(2L+1)'
	      do l2=0,lmpb(isort)
	        write(iun,'(i2,5x,i3,8x,i4)')l2,ntle_pb(l2,isort),
     &		                             (l2+l2+1)*ntle_pb(l2,isort)
	      enddo
	    enddo
	  endif
	endif
	deallocate(s0,evec,eval,pp,pp1,qq,qq1,pb,work)
      if(key==1) call timel('**** prod_b2 finished **************')
      end
