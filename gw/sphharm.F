*****************************************************************
*                                                               *
*   package of routines for evaluating spherical harmonics      *
*                                                               *
*            y(l,m)=bett(l,m)*yy(l,m)                           *
*                                                               *
*****************************************************************
*                                                               *
*                        / cos(m*fi) ,m>0                       *
*      yy(l,m)=plj(l,n)*<                  ; n=abs(m)           *
*                        \ sin(n*fi) ,m<0                       *
*                                                               *
*      plj(l,n)-augemented legandre polinoms                    *
*                                                               *
*     parameters :                                              *
*       x,y,z - components of the unit vector                   *
*       ils - maximum value of the orbit number                 *
*       ylm - resulting array                                   *
*                                                               *
*  obtained from K.Hummler in Max-Planck-Institut, Germany      *
*****************************************************************
      subroutine sphharm(x0,y0,z0,ils,ylm)
	use manager_mod
	implicit none
	integer, intent(in) :: ils
	real*8, intent(in) :: x0,y0,z0
      real*8, intent(out) :: ylm((ils+1)**2)
	integer n,nt,l,lp1,lm1,lav,lavpl,lavml,lmm,m,lavpm,lavmm,mp1,l0,lm
	real*8 x,y,z,st,z2,r2
      real*8, allocatable :: yy(:),c(:),s(:),p(:,:)
      allocate(yy((maxylm+1)**2),c(maxylm+1),s(maxylm+1),
     &	       p(maxylm+1,maxylm+1))
      c(1)=1.d0
      s(1)=0.d0
      p(1,1)=1.d0
      p(2,2)=1.d0
      n = (ils+1)**2
      yy(1) = 1.d0
      x = x0
      y = y0
      z = z0
      c(2)=x
      s(2)=y
      p(2,1)=z
      st = x*x + y*y
      z2 = z*z
      r2 = st+z2
      if (n>=2) then
        yy(2) = y
        yy(3) = z
        yy(4) = x
        nt = 1
        do  l = 2, ils
          lp1 = l+1
          lm1 = l-1
          lav = l*lp1 + 1
          p(lp1,1) = ((l+lm1)*z*p(l,1) - lm1*r2*p(lm1,1)) / l
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
              p(lp1,mp1)=(r2*(l+m)*p(l,mp1)-st*p(lp1,mp1+1))/(z*(l-m))
              yy(lavpm) = p(lp1,mp1)*c(mp1)
              yy(lavmm) = p(lp1,mp1)*s(mp1)
	      enddo
          else
            do m = 1, lm1
              mp1 = m+1
              lavpm = lav+m
              lavmm = lav-m
              p(lp1,mp1)=((lm1+m)*r2*p(l,m)-(lp1-m)*z*p(lp1,m))/st
              yy(lavpm) = p(lp1,mp1)*c(mp1)
              yy(lavmm) = p(lp1,mp1)*s(mp1)
	      enddo
          endif
	  enddo
      endif
c  inserted by s.s -->
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
      
      
      subroutine sph_harm_cmp(xx,yy,zz,lmax,ylm)
	use atom_mod
	use etot_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
      integer, intent(in) :: lmax
      real*8, intent(in) :: xx,yy,zz
      complex*16, intent(out) :: ylm((lmax+1)**2)
      integer :: l,m,lm,lim,ll,lm1
      real*8 :: sq
      real*8, allocatable :: ylmr(:)
      lim=(lmax+1)**2
      sq=1.d0/sqrt(2.d0)
      allocate(ylmr(lim))
      call sphharm(xx,yy,zz,lmax,ylmr)
      do l=0,lmax
        ll=l*(l+1)+1
        do m=-l,-1
          lm=ll+m
          lm1=ll-m
          ylm(lm)=sq*dcmplx(ylmr(lm1),-ylmr(lm))
        enddo
        ylm(ll)=ylmr(ll)
        do m=1,l
          lm=ll+m
          lm1=ll-m
          ylm(lm)=sq*(-1)**m*dcmplx(ylmr(lm),ylmr(lm1))
        enddo
      enddo   !! over l
      deallocate(ylmr)
      end
      
      
      subroutine bispinor(xx,yy,zz,lmax,bisp)
	use atom_mod
	use etot_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
      integer, intent(in) :: lmax
      real*8, intent(in) :: xx,yy,zz
      complex*16, intent(out) :: bisp(2*(lmax+1)**2,2)
      integer :: l,lm,lim,ll,ist,l2,l1,ii,li,jj,mj,limj,ind
      real*8 :: sq
      complex*16, allocatable :: ylm(:)
      lim=(lmax+1)**2
      sq=1.d0/sqrt(2.d0)
      allocate(ylm(lim))
      call sph_harm_cmp(xx,yy,zz,lmax,ylm)
      bisp=(0.d0,0.d0)
	do l=0,lmax
	  ist=-1
        if(l==0) ist=1
        ll=l*(l+2)
        l2=l+l
        l1=l*(l+1)+1
	  do ii=ist,1,2
          li=(5+ii)/2+2*(l-1)
          jj=l2+ii
          do mj=-jj,jj,2
            call getlimj(limj,l,ii,mj,li,1)
            ind=ll+(mj+3)/2         !!!  IV-37-7
            if(mj+1<=l2) then
              lm=l1+(mj+1)/2
              bisp(limj,1)=-ckg2(ii,1,ind)*ylm(lm)
            endif
            if(mj-1>=-l2) then
              lm=l1+(mj-1)/2
              bisp(limj,2)=ckg2(ii,2,ind)*ylm(lm)
            endif
          enddo
        enddo
      enddo   !! over l
      deallocate(ylm)
      end