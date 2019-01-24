      subroutine rad_eq3(e,z,l,i,nodes,v,q,r,nrad,pnew,qnew,fi,dfi,dl,
     &                   dr,h,p0,q0,p1,work,f_ex,f_exq,isort,ispin,li,
     &                   node)
c	The outward solution for given e

c	Input:
c		pnew = P_old(r)  - for key=1 or 2 this is the first or
c		qnew = Q_old(r)    second derivatives correspondingly
c	  p0,q0 = P,Q for key=1,2; work array for key=0
c	  p1,q1 = P_dot,Q_dot for key=2; aren't needed for key < 2
c		nodes = number of nodes
c	Output:
c		pnew = P_new(r)
c		qnew = Q_new(r)
c		fi =  P(S)/S
c		dfi = (P(S)/S)'
c		nod = number of nodes
c		dl = log derivative
c		q = <P!P>
	use manager_mod
	use units_mod
      implicit none
	integer, intent(in) :: l,i,nrad,isort,ispin,li,nodes
	real*8, intent(in) :: z,v(0:nrad),r(0:nrad),h,p1(0:nrad),
     &	                    dr(0:nrad)
      integer, intent(out) :: node
      real*8, intent(out) :: q,fi,dfi,dl
      real*8, intent(inout) :: work(0:nrad),p0(0:nrad),q0(0:nrad),
     &                         f_ex(0:nrad),f_exq(0:nrad),e,
     &                         pnew(0:nrad),qnew(0:nrad)
	integer :: kap,it
	real*8 :: c2,b,akap,fac,dp,wk,alf,dfi_big,fff,e0,e2,e1,de
	real*8, allocatable :: pn(:),qn(:),fex(:),fexq(:)
	allocate(pn(0:nrad),qn(0:nrad))
	allocate(fex(0:nrad),fexq(0:nrad))
      c2=clight**2
      b=r(nrad)/(dexp(h*nrad-h)-1.d0)
	if(irel.eq.0) then
	  akap=dfloat(l+1)
      else if(irel.eq.1) then
        akap=dfloat(l*(l+1))
      else if(irel.eq.2) then
        kap=0
        if(i.lt.0) kap=l
        if(i.gt.0) kap=-l-1
        akap=dfloat(kap)
	endif
	it=1
	fff=1.d0
	do while(it<=10.and.fff>1.d-12)
	  fex=f_ex
	  fexq=f_exq
        call RADSCH(e,z,h,b,l,nrad,r,dr,pnew,qnew,f_exq,f_ex,v,1,akap,
     &              1,pnew(1),qnew(1),p0,q0,0,nrad,node,q,work,p1,irel)
        if(node.ne.nodes) then
          e0=e
          if(node.gt.nodes) e2=e
          if(node.lt.nodes) e1=e
          if(node.gt.nodes) e=min(e*1.05d0,e-0.03d0)
          if(node.lt.nodes) e=max(e*0.95d0,e+0.03d0)
          de=e-e0
          goto 1
        endif
c -----  normalize Pnew, Qnew -------
        fac=1.d0/dsqrt(q)
	  pnew=fac*pnew
	  qnew=fac*qnew
	  call f_ex_new(f_ex,f_exq,isort,pnew,qnew,ispin,li,nrad)
	  fff=maxval(abs(f_ex-fex))
	  f_ex=0.5d0*f_ex+0.5d0*fex
	  f_exq=0.5d0*f_exq+0.5d0*fexq
        it=it+1
1       continue        
	enddo
	fi=pnew(nrad)/r(nrad)
	if(irel.eq.0) then
	  dp=(l+1)/r(nrad)*pnew(nrad)+qnew(nrad)
	  dfi=dp/r(nrad)-pnew(nrad)/r(nrad)**2
	  dl=dfi/fi*r(nrad)
	else if(irel.eq.1) then
        wk=v(nrad)-(z+z)/r(nrad)-e
        alf=1.d0-wk/c2
	  fac=dfloat(l)*(l+1.d0)/alf/r(nrad)**2
	  dp=pnew(nrad)/r(nrad)+alf*qnew(nrad)
	  dp=dp-f_exq(nrad)/c2
	  dfi=dp/r(nrad)-pnew(nrad)/r(nrad)**2
	  dl=dfi/fi*r(nrad)
	else if(irel.eq.2) then
        wk=v(nrad)-(z+z)/r(nrad)-e
        alf=1.d0-wk/c2
	  dp=-akap*pnew(nrad)/r(nrad)+alf*qnew(nrad)
	  dp=dp-f_exq(nrad)/c2
	  dfi_big=dp/r(nrad)-pnew(nrad)/r(nrad)**2
	  dl=dfi_big/fi*r(nrad)
	  if(rel_interst) dfi=qnew(nrad)/r(nrad)
	  if(.not.rel_interst) dfi=dfi_big
	endif
	deallocate(pn,qn,fex,fexq)
      end
