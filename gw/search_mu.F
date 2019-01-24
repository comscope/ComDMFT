      subroutine search_mu(x1,y0,x0)
c ------ Search for Chemical Potential --------------------------
c Adapted from Yu.A.Uspensky-s code
c     key=0 for free ( or LDA) electrons
c     key=1 for exchange only
c     key=2 for full interaction
	use manager_mod
	use parallel_mod
      implicit none
      logical left,right
      integer :: key,itmax,it
      real*8 :: f1,f2,c,d,y,y0,deltx,x1,x0,x,yy,dy,
     &          xl,yyl,dyl,xr,yyr,dyr,eps,h,h1,a,b,dys,alfa,beta
      common/ltrt/x,yy,dy,xl,yyl,dyl,xr,yyr,dyr,right,left
c
      f1(a,b)=2.d0*a/(1.d0+sqrt(1.+4.d0*a*(b-1)))
      f2(c,d)=2.d0/(c+sqrt(c*c+4.d0*(d-c)))
c
c               initiation
	if(ubi=='dft') key=0
	if(ubi==' hf') key=1
	if(ubi==' gw') key=2
	if(ubi=='gwg') key=2
	h=0.2d0
	eps=1.d-9
	itmax=50
      it=1
      h1=1.d0/h
      x=x1
      left=.false.
      right=.false.
c
c       the beginning of the iteration loop.
c
1     call tos(x,y,dy,key)
      yy=y-y0
      call sltrt
      if(.not.left.or..not.right) then
c       the newton-s method.
        if(dy.le.abs(yy*h1)) x=x-sign(h,yy)  !! with limitation
        if(dy.gt.abs(yy*h1)) x=x-yy/dy       !! without limitation
      else if(left.and.right) then
c       the parabola method.
        deltx=xr-xl
        if(deltx.eq.0.)goto 2
        dys=(yyr-yyl)/deltx
        if(abs(yyl).le.yyr) then
          alfa=abs(yyl/deltx)
          if(alfa.le.dyl) then
            a=alfa/dyl
            b=dys/dyl
            x=xl+deltx*f1(a,b)
          else
            c=dyl/alfa
            d=dys/alfa
            x=xl+deltx*f2(c,d)
          end if
        else
          beta=yyr/deltx
          if(beta.le.dyr) then
            a=beta/dyr
            b=dys/dyr
            x=xr-deltx*f1(a,b)
          else
            c=dyr/beta
            d=dys/beta
            x=xr-deltx*f2(c,d)
          end if
        end if
      end if
      if(it.ge.itmax)then
       if(maswrk) write(iun,*)
     &	        ' The exceeding of itmax in chem_pot_search'
         goto 2
      endif
      if(abs(yy).le.eps) go to 2
c
      it=it+1
      go to 1
c
c               the end of the loop.
c
2     if(it.ne.1)then
        x0=x
      else
        x0=x1
      endif
      end
