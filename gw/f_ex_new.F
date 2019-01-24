      subroutine f_ex_new(f_ex,f_exq,isort,pp,qq,ispin,li,nw)
	use atom_mod
	use manager_mod
	use solid_mod
	use units_mod
      implicit none
	integer, intent(in) :: isort,ispin,li,nw
	real*8, intent(in) :: pp(0:nw),qq(0:nw)
	real*8, intent(out) :: f_ex(0:nw),f_exq(0:nw)
	integer :: li1,lf1,l1,lf2,li2,l2,l,ic1,mt1,ir,mt2,ll2,isp
	real*8 :: wh1
	real*8, allocatable :: pval(:),qval(:),pp1(:),f(:),g(:),qq1(:),
     &	                     tmp(:),work(:)
      isp=min(ispin,nspin_1)
	allocate(tmp(0:nw),pp1(0:nw),f(0:nw),g(0:nw),qval(0:nw),qq1(0:nw),
     &	       work(0:nw),pval(0:nw))
	f_ex(0:nw)=0.d0
	f_exq(0:nw)=0.d0
	if(irel_core.le.1) then
	  l=li-1
	else
	  l=li/2
	endif
c ------ core contribution --------------
	do ic1=1,ncor(isort)
	  li1=li_cor(ic1,isort)
	  if(irel_core.le.1) then
	    l1=li1-1
	  else
	    l1=li1/2
	  endif
	  mt1=indcor(ic1,ispin,isort)
	  call dcopy(nrad(isort)+1,pcor_old(mt1),1,pp1(0),1)
	  if(irel.ge.1) call dcopy(nrad(isort)+1,qcor_old(mt1),1,qq1(0),1)
	  do l2=0,min(l+l1,max_rl)
	    call var_int(dr(0,isort),pp(0),qq(0),pp1(0),qq1(0),f(0),g(0),
     &	               tmp(0),r(0,isort),rl(0,l2,isort),work(0),
     &                 nrad(isort)+1,h(isort))
	    wh1=t_x(li,li1,l2)
	    call dcopy(nrad(isort)+1,tmp(0),1,f(0),1)
	    call dscal(nrad(isort)+1,wh1,f(0),1)
	    do ir=0,nrad(isort)
	      f_ex(ir)=f_ex(ir)+f(ir)*pp1(ir)
	    enddo
	    if(irel_core.ge.1) then
	      do ir=0,nrad(isort)
	        f_exq(ir)=f_exq(ir)+f(ir)*qq1(ir)
	      enddo
	    endif
	  enddo   !! over l2
	enddo  !! over ic1
c ------ Valence contribution ----------------------------------------
	do lf2=1,lfun(isort)
	  li2=li_val(lf2,isort)
	  if(irel.le.1) then
	    ll2=li2-1
	  else
	    ll2=li2/2
	  endif
	  mt2=ind_wf(lf2,isort)
	  do ir=0,nrad(isort)
	    pval(ir)=gfun_old(mt2+ir,isp)*r(ir,isort)
	  enddo
	  if(irel.ge.1) then
	    do ir=0,nrad(isort)
	      qval(ir)=gfund_old(mt2+ir,isp)*r(ir,isort)
	    enddo
	  endif
	  do l2=0,min(l+ll2,max_rl)
	    call var_int(dr(0,isort),pval(0),qval(0),pp(0),qq(0),
     &			       f(0),g(0),tmp(0),r(0,isort),rl(0,l2,isort),
     &                 work(0),nrad(isort)+1,h(isort))
	    do lf1=1,lfun(isort)
	      mt1=ind_wf(lf1,isort)
	      do ir=0,nrad(isort)
	        f(ir)=gfun_old(mt1+ir,isp)*r(ir,isort)
	        if(irel.ge.1) g(ir)=gfund_old(mt1+ir,isp)*r(ir,isort)
	      enddo
	      do ir=0,nrad(isort)
	        f_ex(ir)=f_ex(ir)+tmp(ir)*t1_x(l2,lf2,lf1,li)*f(ir)
	      enddo
	      if(irel>=1) then
	        do ir=0,nrad(isort)
	          f_exq(ir)=f_exq(ir)+tmp(ir)*t1_x(l2,lf2,lf1,li)*g(ir)
 	        enddo
	      endif
	    enddo   !! over lf1
	  enddo   !! over l2
	enddo  !! over lf2
	deallocate(tmp,pp1,f,g,qval,qq1,work,pval)
      end