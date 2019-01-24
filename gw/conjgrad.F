      subroutine conjgrad(func,n,x,f,eps,nit)
      implicit none
	external func
      integer, intent(in) :: n,nit
	real*8, intent(in) :: eps
      real*8, intent(inout) :: x(n)
	real*8, intent(out) :: f
	real*8, allocatable :: grad(:),grad0(:),h0(:),h(:),x0(:),xmi(:)
	integer :: iter,mline,i
	real*8 :: stepalf,anorm,anormp,dspmax,anorm0,alf,alfl,alfr,dalfl,
     &	        dalfr,ddot,bet,funcl,difalf,funcr,xx,cubint,dalf,eps1,
     &          fmi
	allocate(grad0(n),grad(n),x0(n),xmi(n))
	allocate(h0(n),h(n))
      anorm0=0.0d0
      stepalf=0.0d0
	dspmax=1.d-2
	x0=x
	mline=0
	fmi=1.d5
	do iter=1,nit+100
	  call func(n,x,f,grad,1)
	  if(f<fmi) then
	    fmi=f
	    xmi(1:n)=x(1:n)
	  endif
        if(mline.eq.0) then
          anorm=ddot(n,grad,1,grad,1)
          if(sqrt(anorm).le.1.d-15) exit
	    if(iter>nit) exit
          stepalf=0.0d0
          if(iter.eq.1) then
	      h=grad
            anormp=0.d0
            do i=1,n
               anormp=max(anormp,abs(grad(i)))
            enddo
            stepalf=dspmax/anormp
            dalfl=-anorm
          else if(iter.gt.1) then
            bet=-anorm/anorm0
	      if(abs(bet)>10.d0) bet=bet/abs(bet)
	      h=grad-bet*h
            dalfl=-ddot(n,grad,1,h,1)
            anormp=0.d0
            do i=1,n
               anormp=max(anormp,abs(grad(i)))
            enddo
            stepalf=min(1.d0,dspmax/anormp)
          endif
          call dcopy(n,grad,1,grad0,1)
          anorm0=anorm
          alfl=0.d0
          funcl=f
          mline=1
          alf=stepalf
          difalf=alf
	    x=x0-alf*h
        else if(mline.eq.1) then
          dalfr=-ddot(n,grad,1,h,1)
          alfr=alf
          funcr=f
          if(dalfr.lt.0.d0) then
            xx=dalfr*(alfr-alfl)/(dalfl-dalfr)
            dalfl=dalfr
            alfl=alf
            funcl=f
            if(xx.lt.0.d0) xx=2.d0*stepalf
            stepalf=min(2.d0*stepalf,xx)
            alf=alf+stepalf
            difalf=stepalf
          else
            eps1=1.d-7*0.5d0*(alfl+alfr)
            alf=cubint(alfl,alfr,funcl,funcr,dalfl,dalfr,eps1)
            mline=2
            difalf=alfr-alf
            if(difalf.lt.0.d0.or.alf<alfl.or.alf>alfr) then
              alf=0.5d0*(alfl+alfr)
              difalf=alf-alfl
            endif
          endif
	    x=x0-alf*h
        else if(mline.eq.2) then
          dalf=-ddot(n,grad,1,h,1)
          if(dalf.lt.0.d0) then
            dalfl=dalf
            alfl=alf
            funcl=f
            eps1=1.d-7*0.5d0*(alfl+alfr)
            alf=cubint(alfl,alfr,funcl,funcr,dalfl,dalfr,eps1)
            difalf=alf-alfl
            if(difalf.lt.0.d0.or.alf<alfl.or.alf>alfr) then
              alf=0.5d0*(alfl+alfr)
              difalf=alf-alfl
            endif
          else
            dalfr=dalf
            alfr=alf
            funcr=f
            eps1=1.d-7*0.5d0*(alfl+alfr)
            alf=cubint(alfl,alfr,funcl,funcr,dalfl,dalfr,eps1)
            difalf=alfr-alf
            if(difalf.lt.0.d0.or.alf<alfl.or.alf>alfr) then
              alf=0.5d0*(alfl+alfr)
              difalf=alf-alfl
            endif
          endif
	    x=x0-alf*h
          if(anormp*difalf.le.eps) then
		  x0=x
		  mline=0
	    endif
        endif
	enddo
	f=fmi
	x(1:n)=xmi(1:n)
      deallocate(grad0,grad,h0,h,x0,xmi)
      end
