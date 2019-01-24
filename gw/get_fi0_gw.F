	subroutine get_fi0_gw
	use atom_mod
	use manager_mod
	use solid_mod
	use units_mod
      implicit none
	integer :: ispin,isort,le1,mt1,mt,le,ir,i,mtpb,jspin
	real*8 :: c2,dqdall
	real*8, allocatable :: pp(:),pp1(:),qq(:),qq1(:),work(:)
	allocate(pp(0:maxnrad),pp1(0:maxnrad),qq(0:maxnrad),
     &	       qq1(0:maxnrad),work(0:maxnrad))
	c2=clight**2
	do jspin=1,nspin
	  do ispin=1,nspin
	    do isort=1,nsort
	      do le1=1,lfun(isort)
	        mt1=ind_wf(le1,isort)
	        do ir=0,nrad(isort)
	          pp1(ir)=gfun(mt1+ir,jspin)*r(ir,isort)
	          if(irel>=1) qq1(ir)=gfund(mt1+ir,jspin)*r(ir,isort)
	        enddo
	        do le=1,lfun(isort)
	          mt=ind_wf(le,isort)
	          do ir=0,nrad(isort)
	            pp(ir)=gfun(mt+ir,ispin)*r(ir,isort)
	            if(irel>=1) qq(ir)=gfund(mt+ir,ispin)*r(ir,isort)
	          enddo
	          do i=1,lfun_pb(isort)
	            mtpb=ind_pb(i,isort)
	            do ir=0,nrad(isort)
	              work(ir)=pp(ir)*pp1(ir)
	              if(irel>=1) work(ir)=work(ir)+qq(ir)*qq1(ir)/c2
                    work(ir)=work(ir)*fun_pb(ir+mtpb)*dr(ir,isort)
                  enddo
                  fi0(le,le1,i,isort,ispin,jspin)=dqdall(h(isort),work,
     &		                                             nrad(isort))
	          enddo   !! over i
	        enddo   !! over le
	      enddo   !! over le1
	    enddo  !! over isort
	  enddo  !! over ispin
	enddo  !! over jspin
	deallocate(pp,pp1,qq,qq1,work)
      end
      
      
      
      subroutine dd_xx_atom(dd,xx,limax)
	use atom_mod
	use manager_mod
	use solid_mod
	use units_mod
      implicit none
      integer, intent(in) :: limax
      real*8, intent(out) :: dd(limax,limax,0:2*maxb),
     &                       xx(limax,limax,0:2*maxb) 
	integer :: i,l,lm1,li1,lget,l1,i1,mj1,lm,li,mj,ii,iii,indx,icg1,
     &	         icg2,icg,lm2,l2
	real*8 :: pi8,div1,div
	complex*16, allocatable :: work(:,:)
	pi8=8.d0*pi
	xx=0.d0
	do lm1=1,nrel*(maxb+1)**2
	  if(irel.le.1) then
	    li1=lget(lm1)+1
	  else if(irel.eq.2) then
	    call getlimj(lm1,l1,i1,mj1,li1,0)
	  endif
	  do lm=1,nrel*(maxb+1)**2
	    if(irel.le.1) then
	      li=lget(lm)+1
	    else if(irel.eq.2) then
	      call getlimj(lm,l,i,mj,li,0)
	    endif
	    if(irel.le.1) then
		    ii = max0(lm1,lm)
		    iii = min0(lm1,lm)
            indx = (ii*(ii-1))/2 + iii
            icg1 = indxcg(indx)
            icg2 = indxcg(indx+1) - 1 
	    else
		    icg1=indxcgr(lm1,lm,1,0)
            icg2=indxcgr(lm1,lm,2,0)
	    endif
          do icg = icg1, icg2
            if(irel.le.1) then
		      lm2 = jcg(icg)
	      else
	        lm2=jcgr(icg)
	      endif
	      l2=lget(lm2)
	      if(irel<=1) xx(li,li1,l2)=xx(li,li1,l2)+cg(icg)**2
	      if(irel==2) xx(li,li1,l2)=xx(li,li1,l2)+abs(cgr(icg))**2
	    enddo
	  enddo
	enddo
c ---------------------------------------------------------
      allocate(work((2*maxb+1)**2,limax))
	work=(0.d0,0.d0)
	do lm=1,nrel*(maxb+1)**2
	  if(irel.le.1) then
	    li=lget(lm)+1
	  else if(irel.eq.2) then
	    call getlimj(lm,l,i,mj,li,0)
	  endif
	  if(irel.le.1) then
	    ii = lm
	    iii = lm
          indx = (ii*(ii-1))/2 + iii
          icg1 = indxcg(indx)
          icg2 = indxcg(indx+1) - 1 
	  else
	    icg1=indxcgr(lm,lm,1,0)
          icg2=indxcgr(lm,lm,2,0)
	  endif
        do icg = icg1, icg2
          if(irel.le.1) then
		    lm2 = jcg(icg)
	    else
	      lm2=jcgr(icg)
	    endif
	    if(irel<=1) work(lm2,li)=work(lm2,li)+cg(icg)
	    if(irel==2) work(lm2,li)=work(lm2,li)+cgr(icg)
	  enddo
	enddo
	dd=0.d0
	do li1=1,limax
	  do li=1,limax
	    do lm2=1,(2*maxb+1)**2
	      l2=lget(lm2)
	      dd(li,li1,l2)=dd(li,li1,l2)+work(lm2,li)*work(lm2,li1)
	    enddo
	  enddo
	enddo
	deallocate(work)
	do li1=1,limax
        div1=0.0d0
	  if(irel.le.1) then
	    div1=dfloat(li1+li1-1)
	  else if(irel.eq.2) then
	    div1=dfloat(((li1+1)/2)*2)
	  endif
	  do li=1,limax
          div=0.0d0
	    if(irel.le.1) then
	      div=dfloat(li+li-1)
	    else if(irel.eq.2) then
	      div=dfloat(((li+1)/2)*2)
	    endif
	    do l2=0,2*maxb
	      dd(li,li1,l2)=dd(li,li1,l2)*pi8/(l2+l2+1)/div/div1
	      xx(li,li1,l2)=xx(li,li1,l2)*pi8/(l2+l2+1)/div/div1
	    enddo
	  enddo
	enddo
      end
