      subroutine step_descend(func,n,x,f,eps,nit)
      implicit none
	external func
      integer, intent(in) :: n,nit
	real*8, intent(in) :: eps
      real*8, intent(inout) :: x(n)
	real*8, intent(out) :: f
	integer :: iter,k,niter0,iter0
	real*8 :: difmax,e0,e1,e2,f0,f1,f2,ee,ff,delt,y1,y2,y3,y4,y5,y6
	niter0=150
	do iter=1,nit
	  difmax=0.d0
	  do k=1,n
	    iter0=0
	    e1=x(k)
	    call func(n,x,f1)
	    delt=max(0.0001d0,x(k)/1000)
	    x(k)=x(k)+delt
	    e2=x(k)
	    call func(n,x,f2)
	    x(k)=x(k)-2.d0*delt
	    e0=x(k)
	    call func(n,x,f0)
1	    if(f0<f1.and.f1<f2) then
	      e2=e1
	      e1=e0
	      f2=f1
	      f1=f0
	      delt=2.d0*delt
	      e0=e1-delt
	      x(k)=e0
	      call func(n,x,f0)
	    else if(f0>f1.and.f1>f2) then
	      e0=e1
	      e1=e2
	      f0=f1
	      f1=f2
	      delt=2.d0*delt
	      e2=e1+delt
	      x(k)=e2
	      call func(n,x,f2)
	    else if(f0>f1.and.f1<f2) then
	      y1=(e0*e0-e1*e1)/2
	      y2=e0-e1
	      y3=f0-f1
	      y4=(e1*e1-e2*e2)/2
	      y5=e1-e2
	      y6=f1-f2
	      ee=(y1*y6-y3*y4)/(y2*y6-y3*y5)
	      delt=abs(ee-e1)
	      x(k)=ee
	      call func(n,x,ff)
	      if(e1<ee.and.ee<e2) then
	        if(ff>f1) then
	          e2=ee
	          f2=ff
	        else if(ff<f1) then
	          e0=e1
	          f0=f1
	          e1=ee
	          f1=ff
	        endif
	      else if(e0<ee.and.ee<e1) then
	        if(ff>f1) then
	          e0=ee
	          f0=ff
	        else if(ff<f1) then
	          e2=e1
	          f2=f1
	          e1=ee
	          f1=ff
	        endif
	      endif
	    endif
	    difmax=max(difmax,abs(f1))
	    iter0=iter0+1
	    if(delt>eps.and.iter0<niter0) goto 1
	  enddo   !! over k
	  if(difmax<eps) exit
	enddo  !! over iter
      f=f1
      end