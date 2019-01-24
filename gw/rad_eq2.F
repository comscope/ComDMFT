      subroutine rad_eq2(e,tol,l,i,nod,v0,fifi,pnew,qnew,f_ex,f_exq,
     &                   work,ie,isort,ispin,nrad0,li,key1,node)
c	The outward solution for Nod and Norm=1
c	Input:
c	  key1=1  - Exchange is included
c		Nodes
c	Output:
c		pnew = P_new(r)
c		qnew = Q_new(r)
	use atom_mod
	use manager_mod
	use units_mod
      implicit none
      integer, intent(in) :: l,i,nrad0,ie,isort,ispin,li,key1,nod
      real*8, intent(in) :: tol,v0(0:nrad0),f_ex(0:nrad0),
     &                      f_exq(0:nrad0)
      real*8, intent(inout) :: pnew(0:nrad0),qnew(0:nrad0),work(0:nrad0)
      integer, intent(out) :: node
      real*8, intent(out) :: e,fifi
      logical :: bis,tw
	integer :: it
      real*8 :: fff,e0,e1,dny1,q,q001,q01,e2,emi,ema,de,e00,ee,q1
      real*8, allocatable :: p5(:),q5(:)
	allocate(p5(0:nrad0),q5(0:nrad0))
	it=0
	fff=1.d0
      q001=1.0d15
      e00=e
      e1=e-5.d0
      e2=e+5.d0
      emi=e
      ema=e
      de=1.0d0
      bis=.false.
	tw=.false.
	do while(it<=100.and.fff>tol)
	  emi=min(emi,e)
	  ema=max(ema,e)
        call rad_eq0(e,z(isort),l,i,node,v0,q,r(0,isort),nrad(isort),
     &	             pnew,qnew,fi(ie,li,isort,ispin),
     &               dfi(ie,li,isort,ispin),dny1,dr(0,isort),0,h(isort),
     &               p5,q5,p5,q5,work,f_ex,f_exq,key1)
        if(node.ne.nod) then
          e0=e
          if(node.gt.nod) e2=e
          if(node.lt.nod) e1=e
          if(node.gt.nod) e=min(e*1.2d0,e-0.1d0)
          if(node.lt.nod) e=max(e*0.85d0,e+0.1d0)
          de=e-e0
          goto 1
        endif
	  if(.not.bis) then
	    if(.not.tw) then
	      de=-max(0.03d0,abs(e)/100.d0)
	      e00=e
	      q001=q-1.d0
	      tw=.true.
	    else
	      e0=e
	      q01=q-1.d0
	      if(q01*q001.lt.0.d0) then
	        bis=.true.
	        ee=0.5d0*(e0+e00)
	        de=ee-e
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
	  if(e+de.ge.e2) then
	    e=0.5d0*(e+e2)
	  else if(e+de.le.e1) then
	    e=0.5d0*(e+e1)
	  else
          e=e+de
	  endif
1	  fff=dabs(de/e)
	  if(dabs(e).lt.1.d0) fff=dabs(de)
	  it=it+1
	enddo
	fifi=q
	deallocate(p5,q5)
      end
