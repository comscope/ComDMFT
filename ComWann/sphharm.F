*****************************************************************
*                                                                 
*     package of routines for evaluating real spherical harmonics 
*                                                                 
*     normalization condition:                                    
*     int_d\omega y_{l,m}(omega) y_{l,m}(omega)
*     =delta_{l,l'}delta_{m,m'}                
*     y(l,m)=bett(l,m)*yy(l,m)                 
*                                              
*****************************************************************
*                                                              *
*                       / cos(m*fi) ,m>0                       *
*     yy(l,m)=plj(l,n)*<                  ; n=abs(m)           *
*                       \ sin(n*fi) ,m<0                       *
*                                                              *
*     plj(l,n)-augemented legandre polinoms                    *
*                                                              *
*     parameters :                                             *
*     x,y,z - components of the unit vector                   *
*     ils - maximum value of the orbit number                 *
*     ylm - resulting array                                   *
*                                                             *
*     obtained from K.Hummler in Max-Planck-Institut, Germany      *
*****************************************************************
      subroutine sphharm(x0,y0,z0,ils,ylm)
      use comwann_mod
      implicit none
      integer, intent(in) :: ils
      double precision, intent(in) :: x0,y0,z0
      double precision, intent(out) :: ylm((ils+1)**2)
      integer n,nt,l,lp1,lm1,lav,lavpl,lavml,lmm,m,lavpm,lavmm,mp1,l0,lm
      double precision x,y,zr,st,z2,r2
      double precision, allocatable :: yy(:),c(:),s(:),p(:,:)
      allocate(yy((maxylm+1)**2),c(maxylm+1),s(maxylm+1),
     &  p(maxylm+1,maxylm+1))
      yy=0.0d0
      c=0.0d0
      s=0.0d0
      p=0.0d0

      c(1)=1.d0
      s(1)=0.d0
      p(1,1)=1.d0
      p(2,2)=1.d0
      n = (ils+1)**2
      yy(1) = 1.d0
      x = x0
      y = y0
      zr = z0
      c(2)=x
      s(2)=y
      p(2,1)=zr
      st = x*x + y*y
      z2 = zr*zr
      r2 = st+z2
      if (n>=2) then
        yy(2) = y
        yy(3) = zr
        yy(4) = x
        nt = 1
        do  l = 2, ils
          lp1 = l+1
          lm1 = l-1
          lav = l*lp1 + 1
          p(lp1,1) = ((l+lm1)*zr*p(l,1) - lm1*r2*p(lm1,1)) / l
          yy(lav) = p(lp1,1)
          nt = nt+2
          p(lp1,lp1) = p(l,l)*nt
          c(lp1) = x*c(l) - y*s(l)
          s(lp1) = x*s(l) + y*c(l)
          lavpl = lav+l
          yy(lavpl) = p(lp1,lp1)*c(lp1)
          lavml = lav-l
          yy(lavml) = p(lp1,lp1)*s(lp1)
          if (st<=z2) then
            do lmm = 1, lm1
              m = l-lmm
              lavpm = lav+m
              lavmm = lav-m
              mp1 = m+1
              p(lp1,mp1)=(r2*(l+m)*p(l,mp1)-st*p(lp1,mp1+1))/(zr*(l-m))
              yy(lavpm) = p(lp1,mp1)*c(mp1)
              yy(lavmm) = p(lp1,mp1)*s(mp1)
            enddo
          else
            do m = 1, lm1
              mp1 = m+1
              lavpm = lav+m
              lavmm = lav-m
              p(lp1,mp1)=((lm1+m)*r2*p(l,m)-(lp1-m)*zr*p(lp1,m))/st
              yy(lavpm) = p(lp1,mp1)*c(mp1)
              yy(lavmm) = p(lp1,mp1)*s(mp1)
            enddo
          endif
        enddo
      endif
c     inserted by s.s -->
      do l=0,ils
        l0=l*(l+1)+1
        ylm(l0)=yy(l0)*bett(l0)
        do m=1,l
          lm=l0+m
          lmm=l0-m
          ylm(lm)=yy(lm)*bett(lm)
          ylm(lmm)=yy(lmm)*bett(lmm)
        enddo
      enddo
      deallocate (yy,c,s,p)
      end
