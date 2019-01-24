      subroutine sltrt
c       the classification of x into left and right points.
      implicit none
      logical left,right
      real*8 :: x,yy,dy,xl,yyl,dyl,xr,yyr,dyr
      common/ltrt/x,yy,dy,xl,yyl,dyl,xr,yyr,dyr,right,left
c
      if(yy.le.0.) then
        xl=x
        yyl=yy
        dyl=dy
        left=.true.
      else
        xr=x
        yyr=yy
        dyr=dy
        right=.true.
      end if
      end