      subroutine search_mu_0(x1,y0,x0)
      use manager_mod
      use parallel_mod
      implicit none
      real*8, intent(in) :: x1,y0
      real*8, intent(out) :: x0
      integer :: key,itmax,it,m,i
      real*8 :: y,xx(4),yy(4),dy,eps,h,h1,x,yy1
      logical ok_fint
c
c     Comments:
c
c     1. Subroutine tos only produces values in y (the number of
c        electrons), dy (the gradient) is not calculated for GW or QP,
c        and cannot be relied upon.
c
c     2. The subroutine fint does not always produce the solution you
c        want. If the solution x is outside the domain [xx(m):xx(m+1)]
c        then use bisection instead.
c
c               initiation
      if(ubi=='dft') key=0
      if(ubi==' hf') key=1
      if(ubi==' qp') key=3
      if(ubi==' gw') key=2
      if(maswrk) then
        write(iun,*)' Chemical potential search'
        if(key==0) then
          write(iun,*)' Chem_pot_LDA Tos Dtos:'
        else if(key==1) then
          write(iun,*)' Chem_pot_x Tos Dtos:'
        else if(key==2) then
          write(iun,*)' Chem_pot_GW   Tos_X   Tos_C   Tos_SUM:'
        else if(key==3) then
          write(iun,*)' Chem_pot_QP Tos:'
        endif
      endif
      h=0.2d0
      eps=1.d-12
      itmax=85
      it=1
      h1=1.d0/h
      x=x1
c
c       the beginning of the iteration loop.
c
3     xx(1)=x
      call tos(xx(1),y,dy,key)
      yy(1)=y-y0
      if(abs(yy(1)).le.eps) go to 2
      m=0
      if(yy(1)<0) m=1
      x=xx(1)-sign(0.001d0,yy(1))
      call tos(x,y,dy,key)
      yy1=y-y0
      if(abs(yy1).le.eps) go to 2
      if(abs(yy(1)-yy1)<eps) then
        if(yy1>0.d0) x=x-0.1d0
        if(yy1<0.d0) x=x+0.1d0
        if(it.ge.itmax) then
          if(maswrk) write(iun,*)
     &              ' The exceeding of itmax in chem_pot_search'
         goto 2
        endif
        it=it+1
        goto 3
      endif
      if(m==0) then
        xx(2)=xx(1)
        yy(2)=yy(1)
        xx(1)=x
        yy(1)=yy1
        if(yy(1)<0) m=1
      else
        xx(2)=x
        yy(2)=yy1
        if(yy(2)<0) m=2
      endif
c ------------------------- Two points ---------------------------
      if(yy(1)*yy(2)<0) then
        call fint(xx,yy,2,0.d0,x)
        if(x<xx(m).or.x>xx(m+1)) x=0.5d0*(xx(m)+xx(m+1))
      else
        dy=(yy(2)-yy(1))/(xx(2)-xx(1))
c       the newton-s method.
        if(yy(2)<0) then
          if(dy.le.abs(yy(2)*h1)) x=xx(2)-sign(h,yy(2))  !! with restr.
          if(dy.gt.abs(yy(2)*h1)) x=xx(2)-yy(2)/dy       !! without
        else if(yy(1)>0) then
          if(dy.le.abs(yy(1)*h1)) x=xx(1)-sign(h,yy(1))  !! with restr.
          if(dy.gt.abs(yy(1)*h1)) x=xx(1)-yy(1)/dy       !! without
        endif
      endif
      call tos(x,y,dy,key)
      yy1=y-y0
      if(abs(yy1).le.eps) go to 2
      do i=3,m+2,-1
        xx(i)=xx(i-1)
        yy(i)=yy(i-1)
      enddo
      xx(m+1)=x
      yy(m+1)=yy1
      if(yy1<0) m=m+1
c ------------------------- Three points -----------------------------
      if(yy(1)*yy(3)<0) then
        call fint(xx,yy,3,0.d0,x)
        if(x<xx(m).or.x>xx(m+1)) x=0.5d0*(xx(m)+xx(m+1))
      else
c       the newton-s method.
        if(yy(3)<0) then
          dy=(yy(3)-yy(2))/(xx(3)-xx(2))
          if(dy.le.abs(yy(3)*h1)) x=xx(3)-sign(h,yy(3))  !! with restr.
          if(dy.gt.abs(yy(3)*h1)) x=xx(3)-yy(3)/dy       !! without
        else if(yy(1)>0) then
          dy=(yy(2)-yy(1))/(xx(2)-xx(1))
          if(dy.le.abs(yy(1)*h1)) x=xx(1)-sign(h,yy(1))  !! with restr.
          if(dy.gt.abs(yy(1)*h1)) x=xx(1)-yy(1)/dy       !! without
        endif
      endif
      call tos(x,y,dy,key)
      yy1=y-y0
      if(abs(yy1).le.eps) go to 2
      do i=4,m+2,-1
        xx(i)=xx(i-1)
        yy(i)=yy(i-1)
      enddo
      xx(m+1)=x
      yy(m+1)=yy1
      if(yy1<0) m=m+1
c ------------------------- Four points --------------------------
1     if(yy(1)*yy(4)<0) then
        if (ok_fint(xx,yy,4,0.0d0)) then
          call fint(xx,yy,4,0.d0,x)
          if(x<xx(m).or.x>xx(m+1)) x=0.5d0*(xx(m)+xx(m+1))
        else
          x = 0.5d0*(xx(m)+xx(m+1))
        endif
      else
c       the newton-s method.
        if(yy(4)<0) then
          dy=(yy(4)-yy(3))/(xx(4)-xx(3))
          if(dy.le.abs(yy(4)*h1)) x=xx(4)-sign(h,yy(4))  !! with restr.
          if(dy.gt.abs(yy(4)*h1)) x=xx(4)-yy(4)/dy       !! without
        else if(yy(1)>0) then
          dy=(yy(2)-yy(1))/(xx(2)-xx(1))
          if(dy.le.abs(yy(1)*h1)) x=xx(1)-sign(h,yy(1))  !! with restr.
          if(dy.gt.abs(yy(1)*h1)) x=xx(1)-yy(1)/dy       !! without
        endif
      endif
      call tos(x,y,dy,key)
      yy1=y-y0
      if(m<2) then
        do i=4,m+2,-1
          xx(i)=xx(i-1)
          yy(i)=yy(i-1)
        enddo
        xx(m+1)=x
        yy(m+1)=yy1
      else if(m==2) then
        if(yy1<0.d0) then
          xx(1)=xx(2)
          yy(1)=yy(2)
          xx(2)=x
          yy(2)=yy1
        else if(yy1>0.d0) then
          xx(4)=xx(3)
          yy(4)=yy(3)
          xx(3)=x
          yy(3)=yy1
        endif
      else if(m>2) then
        do i=1,m-1
          xx(i)=xx(i+1)
          yy(i)=yy(i+1)
        enddo
        xx(m)=x
        yy(m)=yy1
      endif
      m=count(mask=yy<0)
      if(it.ge.itmax)then
         if(maswrk) write(iun,*)
     &              ' The exceeding of itmax in chem_pot_search'
         goto 2
      endif
      if(abs(yy1).le.eps) go to 2
c
      it=it+1
      go to 1
c
c               the end of the loop.
c
2     x0=x
      end
