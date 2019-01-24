      subroutine rad_eq0(e,z,l,i,nod,v,q,r,nrad,pnew,qnew,fi,dfi,dl,dr,
     &                   key,h,p0,q0,p1,q1,work,f_ex,f_exq,key2)
c	The outward solution for given e
c	key=0  - for function
c	key=1  - for first energy derivative
c	key=2  - for second energy derivative
c	key2=1 - Exchange is included

c	Input:
c		pnew = P_old(r)  - for key=1 or 2 this is the first or
c		qnew = Q_old(r)    second derivatives correspondingly
c	  p0,q0 = P,Q for key=1,2; work array for key=0
c	  p1,q1 = P_dot,Q_dot for key=2; aren't needed for key < 2
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
	integer, intent(in) :: l,i,nrad,key,key2
	real*8, intent(in) :: z,v(0:nrad),r(0:nrad),h,p1(0:nrad),
     &	                    q1(0:nrad),dr(0:nrad),e,f_ex(0:nrad),
     &                      f_exq(0:nrad)
      real*8, intent(inout) :: work(0:nrad),p0(0:nrad),q0(0:nrad)
      integer, intent(out) :: nod
      real*8, intent(out) :: pnew(0:nrad),qnew(0:nrad),q,fi,dfi,dl
	integer :: ir,kap,key1
	real*8 :: c2,b,akap,fac,dp,wk,alf,dfi_big,dqdall
	real*8, allocatable :: pn(:),qn(:)
	allocate(pn(0:nrad),qn(0:nrad))
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
	if(key2==0) key1=0
	if(key2==1) key1=1
	if(key<=1) then
        call RADSCH(e,z,h,b,l,nrad,r,dr,pn,qn,f_exq,f_ex,v,key2,akap,
     &              key1,pnew(1),qnew(1),p0,q0,key,nrad,nod,q,work,p1,
     &              irel)
	else if(key.eq.2) then
        call RADSCH(e,z,h,b,l,nrad,r,dr,pn,qn,f_exq,f_ex,v,key2,akap,
     &              key1,pnew(1),qnew(1),p1,q1,key,nrad,nod,q,work,p1,
     &              irel)
	endif
	if(key2==0) then
	  if(key.eq.0) then
c -----  normalize Pn, Qn -------
          fac=1.d0/dsqrt(q)
	    pnew=fac*pn
	    qnew=fac*qn
	    q=1.d0
        else if(key>0) then
c -----  Require the orthogonality with P0 -------
	    do ir=0,nrad
	      work(ir)=p0(ir)*pn(ir)*dr(ir)
          enddo
	    if(irel.ge.1) then
	      do ir=0,nrad
	        work(ir)=work(ir)+q0(ir)*qn(ir)*dr(ir)/c2
            enddo
	    endif
          fac=dqdall(h,work,nrad)
	    pnew=pn-fac*p0
	    qnew=qn-fac*q0
	    q=q-fac*fac
	  endif
	else
	  pnew=pn
	  qnew=qn
	endif
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
	  if(key.eq.1) dp=dp+q0(nrad)/c2
	  if(key.eq.2) dp=dp+2.d0*q1(nrad)/c2
	  if(key2==1) dp=dp-f_exq(nrad)/c2
	  dfi=dp/r(nrad)-pnew(nrad)/r(nrad)**2
	  dl=dfi/fi*r(nrad)
	else if(irel.eq.2) then
        wk=v(nrad)-(z+z)/r(nrad)-e
        alf=1.d0-wk/c2
	  dp=-akap*pnew(nrad)/r(nrad)+alf*qnew(nrad)
	  if(key.eq.1) dp=dp+q0(nrad)/c2
	  if(key.eq.2) dp=dp+2.d0*q1(nrad)/c2
	  if(key2==1) dp=dp-f_exq(nrad)/c2
	  dfi_big=dp/r(nrad)-pnew(nrad)/r(nrad)**2
	  dl=dfi_big/fi*r(nrad)
	  if(rel_interst) dfi=qnew(nrad)/r(nrad)
	  if(.not.rel_interst) dfi=dfi_big
	endif
	deallocate(pn,qn)
      end
