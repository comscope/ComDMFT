      subroutine rad_eq(eb1,eb2,e,tol,z,l,i,nod,val,slo,v,q,r,nrad,pnew,
     &	                qnew,f_ex,fi00,dfi00,kc,nit,de,nre,f_exq,dr,
     &                  key1,h,dl,work,acc,node,irel0)
c	The Bound solution for given L and NOD
c	key1 = 1  Exchange is included
c	Input:
c		val = P(S)
c		slo = P'(S)
c		nod = number of nodes
c		f_ex =
c		f_exq = 
c	Output:
c		pnew = P_new(r)
c		qnew = Q_new(r)
c		fi =  P(S)/S
c		dfi = (P(S)/S)'
c		dl = log derivative
c		q = <P!P>
	use manager_mod
	use units_mod
      implicit none
      integer, intent(in) :: irel0,nrad,nod,key1,l,i
      real*8, intent(in) :: v(0:nrad),val,slo,h,r(0:nrad),dr(0:nrad),
     &                      f_exq(0:nrad),f_ex(0:nrad),z,tol,eb1,eb2
	integer, intent(out) :: node,kc,nit
      real*8, intent(out) :: acc,dl,de,dfi00,fi00,q
	integer :: nre,nitmax,k,nod2,nod1,ir
      real*8 :: pnew(0:nrad),qnew(0:nrad),work(0:nrad),e,c2,e1,e2,akap,
     &          pn,qn,alf,ph_in,qh_in,ratio,q1,q2,pi_in,qi_in,e00,
     &          q001,q01,ee,fff,b,p0,q0,p1,qq1,e0,fac,dp,dfi_big,h12,
     &          emi,ema
	logical bis,tw
	real*8, allocatable :: g1(:),g2(:)
      nitmax=280
      c2=clight**2
      q001=1.0d15
      de=1.0d0
      nit=0
c 	if(l==0.and.nod==2) nit=1
      e1=eb1
      e2=eb2
      e00=e
      b=r(nrad)/(dexp(h*nrad-h)-1.d0)
c ----- start iterations to find energy ------
      bis=.false.
	tw=.false.
	emi=e
	ema=e  
  1   nit=nit+1
      if(nit.gt.nitmax) return
	if(e.le.e1.or.e.ge.e2) e=.5d0*(e1+e2)
	emi=min(emi,e)
	ema=max(ema,e)
	if(key1==0.or.nre==nrad) then
	  call get_nre(e,nre,l,r,v,z,nrad)
	  nre=min(nre,nrad)
	endif
      call pn_qn(irel0,e,z,l,i,val,slo,v(nre),r(nre),pn,qn,f_exq(nre),
     &           key1,akap)
      kc=nre*2/3
	if(mod(kc,2).eq.1) kc=kc-1
	if(kc.ge.nre-20) kc=nre-20
c -------- Homogeneous equation ---------------------------------------
c	------- Inward integration ----------------------------------
      call RADSCH_b(e,z,h,b,l,kc,r,dr,pnew,qnew,f_exq,f_ex,v,0,
     &	            akap,pn,qn,p0,q0,0,nre,nod2,q2,work,irel0)
	ph_in=pnew(kc)
	qh_in=qnew(kc)
c	Outward integration ----------------------------------------------
      call RADSCH(e,z,h,b,l,kc,r,dr,pnew,qnew,f_exq,f_ex,v,0,
     &	          akap,0,p1,qq1,p0,q0,0,nre,nod1,q1,work,p0,irel0)
	if(key1.eq.0) then  !! id est exchange is not included
	  ratio=pnew(kc)/ph_in
	  do ir=0,kc
	    pnew(ir)=pnew(ir)/ratio
	    qnew(ir)=qnew(ir)/ratio
	  enddo
	  q1=q1/ratio**2
	  q=q1+q2
	  node=nod1+nod2
c	  call rad_hf_check(pnew,qnew,nre,node,nrad)
	else if(key1.eq.1) then  !! id est exchange is included
c -------- Inhomogeneous equation --------------------------------------
	  allocate(g1(0:nre))
	  allocate(g2(0:nre))
c	------- Inward integration ----------------------------------
        call RADSCH_b(e,z,h,b,l,kc,r,dr,g1,g2,f_exq,f_ex,v,1,
     &	              akap,pn,qn,p0,q0,0,nre,nod2,q2,work,irel0)
     	  pi_in=g1(kc)
     	  qi_in=g2(kc)
c	Outward integration ----------------------------------------------
        call RADSCH(e,z,h,b,l,kc,r,dr,g1,g2,f_exq,f_ex,v,1,
     & 	            akap,0,p1,qq1,p0,q0,0,nre,nod1,q1,work,p0,irel0)
	  call inhomo(ph_in,pi_in,qh_in,qi_in,pnew,qnew,g1,g2,
     &	            node,q,nre,kc,dr,h)
	  deallocate(g1,g2)
	  call rad_hf_check(pnew,qnew,nre,node,nrad)
 	endif
 	if(key1==0) then
        if(node.ne.nod) then
          if(node.gt.nod) e2=e
          if(node.lt.nod) e1=e
          if(node.gt.nod) e=min(e*1.5d0,e-0.1d0)
          if(node.lt.nod) e=max(e*0.7d0,e+0.1d0)
          goto 1
        endif
      endif
	if(key1.eq.0) then
	  de=ph_in*(qnew(kc)-qh_in)/q
	else
	  if(.not.bis) then
	    if(.not.tw) then
	      de=-max(0.1d0,abs(e)/2000.d0)
            if(node<nod) de=-de
	      e00=e
	      q001=q-1.d0
	      tw=.true.
	    else
	      e0=e
	      q01=q-1.d0
	      if(q01*q001.lt.0.d0) then
	        if(node==nod) then
	          bis=.true.
	          ee=0.5d0*(e0+e00)
	          de=ee-e
	        else if(node>nod) then
	          de=-max(0.1d0,abs(e)/100.d0)
	          e=emi
	          tw=.false.
	        else if(node<nod) then
	          de=max(0.1d0,abs(e)/100.d0)
	          e=ema
	          tw=.false.
	        endif
	      else if(abs(q01).lt.abs(q001)) then
	        e00=e0
	        q001=q01
	      else
	        e=e-de
	        de=-de/2
            endif
          endif
	  else
	    q1=q-1.d0
	    if(q1*q001.lt.0.d0) then
	      q01=q1
	      e0=e
	      e=0.5d0*(e00+e0)
	      de=e-e0
            e=e0
          else
	      q001=q1
	      e00=e
	      e=0.5d0*(e00+e0)
	      de=e-e00
            e=e00 
          endif
        endif
	endif 
	if(e+de.ge.e2) then
	  e=0.5d0*(e+e2)
	else if(e+de.le.e1) then
	  e=0.5d0*(e+e1)
	else
        e=e+de
	endif
	fff=dabs(de/e)
	if(dabs(e).lt.1.d0) fff=dabs(de)
      if(fff.gt.tol.and.nit.lt.nitmax) goto 1
c -----  normalize g -------
      fac=1.d0/dsqrt(q)
      do k=1,nre
        pnew(k)=pnew(k)*fac
        qnew(k)=qnew(k)*fac
      enddo
      do k=nre+1,nrad
        pnew(k)=0.d0
        qnew(k)=0.d0
      enddo
      if(key1==0) q=1.d0
	fi00=pnew(nrad)/r(nrad)
	h12=1.d0/12.d0/h
	if(irel0.eq.0) then
	  dp=(l+1)/r(nrad)*pnew(nrad)+qnew(nrad)
	  dfi00=dp/r(nrad)-pnew(nrad)/r(nrad)**2
	  if(fi00/=0.d0) dl=dfi00/fi00*r(nrad)
	  if(fi00==0.d0) dl=0.d0
	else if(irel0.eq.1) then
        work(nrad)=v(nrad)-(z+z)/r(nrad)-e
        alf=1.d0-work(nrad)/c2
	  fac=dfloat(l)*(l+1.d0)/alf/r(nrad)**2
	  dp=pnew(nrad)/r(nrad)+alf*qnew(nrad)
	  if(key1.eq.1) dp=dp-f_exq(nrad)/c2
	  dfi00=dp/r(nrad)-pnew(nrad)/r(nrad)**2
	  if(fi00/=0.d0) dl=dfi00/fi00*r(nrad)
	  if(fi00==0.d0) dl=0.d0
	else if(irel0.eq.2) then
        work(nrad)=v(nrad)-(z+z)/r(nrad)-e
        alf=1.d0-work(nrad)/c2
	  dp=-akap*pnew(nrad)/r(nrad)+alf*qnew(nrad)
	  if(key1.eq.1) dp=dp-f_exq(nrad)/c2
	  dfi_big=dp/r(nrad)-pnew(nrad)/r(nrad)**2
	  if(fi00/=0.d0) dl=dfi_big/fi00*r(nrad)
	  if(fi00==0.d0) dl=0.d0
	  dfi00=qnew(nrad)/r(nrad)
	endif
c ------ Check of numerical accuracy ----------------
c ------ by comparing numerical derivatives with the right parts ----
      acc=0.d0
      do ir=1,nrad
        if(ir==0) then
          p0=h12*(-25.d0*pnew(0)+48.d0*pnew(1)-36.d0*pnew(2)
     &            +16.d0*pnew(3)-3.d0*pnew(4))/dr(ir)
          q0=h12*(-25.d0*qnew(0)+48.d0*qnew(1)-36.d0*qnew(2)
     &            +16.d0*qnew(3)-3.d0*qnew(4))/dr(ir)
        else if(ir==1) then
          p0=h12*(-3.d0*pnew(0)-10.d0*pnew(1)+18.d0*pnew(2)
     &            -6.d0*pnew(3)+pnew(4))/dr(ir)
          q0=h12*(-3.d0*qnew(0)-108.d0*qnew(1)+18.d0*qnew(2)
     &            -6.d0*qnew(3)+qnew(4))/dr(ir)
        else if(ir==nrad-1) then
          p0=h12*(-pnew(nrad-4)+6.d0*pnew(nrad-3)-18.d0*pnew(nrad-2)
     &            +10.d0*pnew(nrad-1)+3.d0*pnew(nrad))/dr(ir)
          q0=h12*(-qnew(nrad-4)+6.d0*qnew(nrad-3)-18.d0*qnew(nrad-2)
     &            +10.d0*qnew(nrad-1)+3.d0*qnew(nrad))/dr(ir)
        else if(ir==nrad) then
          p0=h12*(3.d0*pnew(nrad-4)-16.d0*pnew(nrad-3)
     &            +36.d0*pnew(nrad-2)-48.d0*pnew(nrad-1)
     &            +25.d0*pnew(nrad))/dr(ir)
          q0=h12*(3.d0*qnew(nrad-4)-16.d0*qnew(nrad-3)
     &            +36.d0*qnew(nrad-2)-48.d0*qnew(nrad-1)
     &            +25.d0*qnew(nrad))/dr(ir)
        else
          p0=h12*(pnew(ir-2)+8.d0*(pnew(ir+1)-pnew(ir-1))-pnew(ir+2))
     &          /dr(ir)
          q0=h12*(qnew(ir-2)+8.d0*(qnew(ir+1)-qnew(ir-1))-qnew(ir+2))
     &          /dr(ir)
        endif
	  if(irel0.eq.0) then
	    p1=(l+1)/r(ir)*pnew(ir)+qnew(ir)
	    q1=(v(ir)-(z+z)/r(ir)-e)*pnew(ir)-(l+1)/r(ir)*qnew(ir)
	    if(key1.eq.1) q1=q1+f_ex(ir)
	  else if(irel0.eq.1) then
          work(ir)=v(ir)-(z+z)/r(ir)-e
          alf=1.d0-work(ir)/c2
	    p1=pnew(ir)/r(ir)+alf*qnew(ir)
	    if(key1.eq.1) p1=p1-f_exq(ir)/c2
	    fac=dfloat(l)*(l+1.d0)/alf/r(ir)**2
	    q1=-qnew(ir)/r(ir)+(work(ir)+fac)*pnew(ir)
	    if(key1.eq.1) q1=q1+f_ex(ir)
	  else if(irel0.eq.2) then
          work(ir)=v(ir)-(z+z)/r(ir)-e
          alf=1.d0-work(ir)/c2
	    p1=-akap*pnew(ir)/r(ir)+alf*qnew(ir)
	    if(key1.eq.1) p1=p1-f_exq(ir)/c2
	    q1=akap*qnew(ir)/r(ir)+work(ir)*pnew(ir)
	    if(key1.eq.1) q1=q1+f_ex(ir)
	  endif
	  acc=max(acc,abs(p1-p0))
	  acc=max(acc,abs(q1-q0))
	enddo
      end
      
      
      
      subroutine pn_qn(irel0,e,z,l,i,val,slo,v,r,pn,qn,f_exq,key1,akap)
	use manager_mod
	use units_mod
      implicit none
      integer, intent(in) :: irel0,l,i,key1
      real*8, intent(in) :: v,e,z,val,slo,r,f_exq
      real*8, intent(out) :: pn,qn,akap
	integer :: kap
      real*8 :: c2,work,alf
      c2=clight**2
	if(irel0.eq.0) then
	  akap=dfloat(l+1)
	  pn=val           !! =P in V-19-3
        qn=slo-akap*val/r     !! =Q in V-19-3
      else if(irel0.eq.1) then
        akap=dfloat(l*(l+1))
        work=v-(z+z)/r-e
        alf=1.d0-work/c2
	  pn=val           !! =P in V-19-3
        if(key1==0) qn=(slo-val/r)/alf 
        if(key1==1) qn=(slo-val/r+f_exq/c2)/alf !in VI-87-28
      else if(irel0.eq.2) then
        kap=0
        if(i.lt.0) kap=l
        if(i.gt.0) kap=-l-1
        akap=kap
        work=v-(z+z)/r-e
        alf=1.d0-work/c2
	  pn=val           !! =P in V-19-3
        if(key1==0) qn=(slo+akap*val/r)/alf
        if(key1==1) qn=(slo+akap*val/r+f_exq/c2)/alf !87-17
	endif
	end
	
	
	
	
	subroutine rad_eq_est(eb1,eb2,e,z,l,i,nod,val,slo,v,r,nrad,
     & 	                    pnew,qnew,f_ex,f_exq,dr,key1,h,work,irel0,
     &                      nre)
c	The Bound solution for given L and NOD
c	key1 = 1  Exchange is included
c	Input:
c		val = P(S)
c		slo = P'(S)
c		nod = number of nodes
c		f_ex =
c		f_exq = 
c	Output:
c		pnew = P_new(r)
c		qnew = Q_new(r)
c		fi =  P(S)/S
c		dfi = (P(S)/S)'
c		dl = log derivative
c		q = <P!P>
	use manager_mod
	use units_mod
      implicit none
      integer, intent(in) :: irel0,nrad,nod,key1,l,i
      real*8, intent(in) :: v(0:nrad),val,slo,h,r(0:nrad),dr(0:nrad),
     &                      f_exq(0:nrad),f_ex(0:nrad),z
      real*8, intent(inout) :: eb1,eb2
	integer :: nre,nod2,nod1,ir,node0,kc,node
      real*8 :: pnew(0:nrad),qnew(0:nrad),work(0:nrad),e,c2,e1,e2,akap,
     &          pn,qn,ph_in,qh_in,ratio,q1,q2,pi_in,qi_in,b,p0,q0,p1,
     &          qq1,qq0,hh,q,e11,e22,hhh
	real*8, allocatable :: g1(:),g2(:)
      c2=clight**2
      e1=eb1
      e2=eb2
      b=r(nrad)/(dexp(h*nrad-h)-1.d0)
      node0=1000
      qq0=1000.d0
      hh=(e2-e1)/200
1     eb1=1000.d0
      eb2=-1000.d0
      e11=1000.d0
      e22=-1000.d0
	do e=e1,e2,hh
	  if(key1==0.or.nre==nrad) then
	    call get_nre(e,nre,l,r,v,z,nrad)
	    nre=min(nre,nrad)
	  endif
        call pn_qn(irel0,e,z,l,i,val,slo,v(nre),r(nre),pn,qn,f_exq(nre),
     &             key1,akap)
        kc=nre*2/3
	  if(mod(kc,2).eq.1) kc=kc-1
	  if(kc.ge.nre-20) kc=nre-20
c -------- Homogeneous equation ---------------------------------------
c	------- Inward integration ----------------------------------
        call RADSCH_b(e,z,h,b,l,kc,r,dr,pnew,qnew,f_exq,f_ex,v,0,
     &	              akap,pn,qn,p0,q0,0,nre,nod2,q2,work,irel0)
	  ph_in=pnew(kc)
	  qh_in=qnew(kc)
c	Outward integration ----------------------------------------------
        call RADSCH(e,z,h,b,l,kc,r,dr,pnew,qnew,f_exq,f_ex,v,0,
     &	            akap,0,p1,qq1,p0,q0,0,nre,nod1,q1,work,p0,irel0)
	  if(key1.eq.0) then  !! id est exchange is not included
	    ratio=pnew(kc)/ph_in
	    do ir=0,kc
	      pnew(ir)=pnew(ir)/ratio
	      qnew(ir)=qnew(ir)/ratio
	    enddo
	    q1=q1/ratio**2
	    q=q1+q2
	    node=nod1+nod2
c	  call rad_hf_check(pnew,qnew,nre,node,nrad)
	  else if(key1.eq.1) then  !! id est exchange is included
c -------- Inhomogeneous equation --------------------------------------
	    allocate(g1(0:nre))
	    allocate(g2(0:nre))
c	------- Inward integration ----------------------------------
          call RADSCH_b(e,z,h,b,l,kc,r,dr,g1,g2,f_exq,f_ex,v,1,
     &	                akap,pn,qn,p0,q0,0,nre,nod2,q2,work,irel0)
     	    pi_in=g1(kc)
     	    qi_in=g2(kc)
c	Outward integration ----------------------------------------------
          call RADSCH(e,z,h,b,l,kc,r,dr,g1,g2,f_exq,f_ex,v,1,
     & 	              akap,0,p1,qq1,p0,q0,0,nre,nod1,q1,work,p0,irel0)
	    call inhomo(ph_in,pi_in,qh_in,qi_in,pnew,qnew,g1,g2,
     &	              node,q,nre,kc,dr,h)
	    deallocate(g1,g2)
	    call rad_hf_check(pnew,qnew,nre,node,nrad)
 	  endif
c 	  if(l==0.and.nod==2) write(2,*)e,node,q
        if(node0==nod.and.node==nod) then
          if(qq0>=1.d0.and.q<=1.d0) then
            eb1=e-hh
            eb2=e
            return
          endif
          if(qq0<=1.d0.and.q>=1.d0) then
            eb1=e-hh
            eb2=e
            return
          endif
          e11=min(e11,e)
          e22=max(e22,e)
        endif
        node0=node
        qq0=q  
 	enddo
 	if(e11>e22) then
 	  if(node<nod) then
 	    hhh=e2-e1
          e1=e2
          e2=e1+hhh*3
        else
          hhh=e2-e1
          e2=e1
          e1=e2-hhh*3
        endif
        goto 1
 	endif
 	if(eb1>eb2) then
 	  eb1=e11
 	  eb2=e22
 	endif
      end
