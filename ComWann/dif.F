

!     Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
!     distributed under the terms of the GNU Lesser General Public
!     License. See the file COPYING for license details.
      
!     BOP
!     !ROUTINE: fderiv_dble
!     !INTERFACE:
      subroutine fderiv_dble(m,n,x,f,g)
!     !INPUT/OUTPUT PARAMETERS:
!     m : order of derivative (in,integer)
!     n : number of points (in,integer)
!     x : abscissa array (in,real(n))
!     f : function array (in,complex(n))
!     g : (anti-)derivative of f (out,complex(n))
!     !DESCRIPTION:
!     Given function $f$ defined on a set of points $x_i$ then if $m\ge 0$ this
!     routine computes the $m$th derivative of $f$ at each point. If $m<0$ the
!     anti-derivative of $f$ given by
!     $$ g(x_i)=\int_{x_1}^{x_i} f(x)\,dx $$
!     is calculated. If $m=-1$ then an accurate integral is computed by fitting
!     the function to a clamped cubic spline_dble. When $m=-3$ the fast but low
!     accuracy trapezoidal integration method is used. Simpson's integration,
!     which is slower but more accurate than the trapezoidal method, is used if
!     $m=-2$.
!     
!     !REVISION HISTORY:
!     Created May 2002 (JKD)
!     EOP
!     BOC
      implicit none
!     arguments
      integer, intent(in) :: m,n
      double precision, intent(in) :: x(n)
      double precision, intent(in) :: f(n)
      double precision, intent(out) :: g(n)
!     local variables
      integer :: i
      double precision ::  x0,x1,x2,dx
!     automatic arrays
      double precision :: cf(3,n)
c$$$      if (n.le.0) then
c$$$        write(*,*)
c$$$        write(*,'("Error(fderiv_dble): invalid number of points : ",I8)') n
c$$$        write(*,*)
c$$$        stop
c$$$      end if
      select case(m)
      case(-3)
!     low accuracy trapezoidal integration
        g(1)=0.d0
        do i=1,n-1
          g(i+1)=g(i)+0.5d0*(x(i+1)-x(i))*(f(i+1)+f(i))
        end do
        return
      case(-2)
!     medium accuracy Simpson integration
        g(1)=0.d0
        do i=1,n-2
          x0=x(i)
          x1=x(i+1)
          x2=x(i+2)
          g(i+1)=g(i)+(x0-x1)*(f(i+2)*(x0-x1)**2
     $      +f(i+1)*(x2-x0)*(x0+2.d0*x1-3.d0*x2)
     $      +f(i)*(x2-x1)*(2.d0*x0+x1-3.d0*x2))/(6.d0*(x0-x2)*(x1-x2))
        end do
        x0=x(n)
        x1=x(n-1)
        x2=x(n-2)
        g(n)=g(n-1)+(x1-x0)*(f(n-2)*(x1-x0)**2
     $    +f(n)*(x1-x2)*(3.d0*x2-x1-2.d0*x0)
     $    +f(n-1)*(x0-x2)*(3.d0*x2-2.d0*x1-x0))/(6.d0*(x2-x1)*(x2-x0))
        return
      case(0)
        g(:)=f(:)
        return
      case(4:)
        g(:)=0.d0
        return
      end select
!     high accuracy integration/differentiation from spline_dble interpolation
      call spline_dble(n,x,f,cf)
      select case(m)
      case(:-1)
        g(1)=0.d0
        do i=1,n-1
          dx=x(i+1)-x(i)
          g(i+1)=g(i)+(((0.25d0*cf(3,i)*dx
     $      +0.3333333333333333333d0*cf(2,i))*dx
     $      +0.5d0*cf(1,i))*dx+f(i))*dx
        end do
      case(1)
        g(:)=cf(1,:)
      case(2)
        g(:)=2.d0*cf(2,:)
      case(3)
        g(:)=6.d0*cf(3,:)
      end select
      return
      end subroutine
      
      
!     Copyright (C) 2011 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
!     This file is distributed under the terms of the GNU Lesser General Public
!     License. See the file COPYING for license details.
      
!     BOP
!     !ROUTINE: spline_dble
!     !INTERFACE:
      subroutine spline_dble(n,x,f,cf)
!     !INPUT/OUTPUT PARAMETERS:
!     n  : number of points (in,integer)
!     x  : abscissa array (in,real(n))
!     f  : input data array (in,double precision(n))
!     cf : cubic spline_dble coefficients (out,double precision(3,n))
!     !DESCRIPTION:
!     Calculates the coefficients of a cubic spline_dble fitted to input data. In other
!     words, given a set of data points $f_i$ defined at $x_i$, where
!     $i=1\ldots n$, the coefficients $c_j^i$ are determined such that
!     $$ y_i(x)=f_i+c_1^i(x-x_i)+c_2^i(x-x_i)^2+c_3^i(x-x_i)^3, $$
!     is the interpolating function for $x\in[x_i,x_{i+1})$. The coefficients are
!     determined piecewise by fitting a cubic polynomial to adjacent points.
!     
!     !REVISION HISTORY:
!     Created November 2011 (JKD)
!     EOP
!     BOC
      implicit none
!     arguments
      integer, intent(in) :: n
      double precision, intent(in) :: x(n)
      double precision, intent(in) :: f(n)
      double precision, intent(out) :: cf(3,n)
!     local variables
      integer ::  i
      double precision :: x0,x1,x2,x3
      double precision :: y0,y1,y2,y3,c1,c2,c3,t0,t1,t2,t3,t4,t5,t6
      if (n.eq.1) then
        cf(:,1)=0.d0
        return
      end if
      if (n.eq.2) then
        cf(1,1)=(f(2)-f(1))/(x(2)-x(1))
        cf(2:3,1)=0.d0
        cf(1,2)=cf(1,1)
        cf(2:3,2)=0.d0
        return
      end if
      if (n.eq.3) then
        x0=x(1)
        x1=x(2)-x0
        x2=x(3)-x0
        y0=f(1)
        y1=f(2)-y0
        y2=f(3)-y0
        t0=1.d0/(x1*x2*(x2-x1))
        t1=x1*y2
        t2=x2*y1
        c1=t0*(x2*t2-x1*t1)
        c2=t0*(t1-t2)
        cf(1,1)=c1
        cf(2,1)=c2
        cf(3,1)=0.d0
        t3=2.d0*c2
        cf(1,2)=c1+t3*x1
        cf(2,2)=c2
        cf(3,2)=0.d0
        cf(1,3)=c1+t3*x2
        cf(2,3)=c2
        cf(3,3)=0.d0
        return
      end if
      y0=f(1)
      y1=f(2)-y0
      y2=f(3)-y0
      y3=f(4)-y0
      x0=x(1)
      x1=x(2)-x0
      x2=x(3)-x0
      x3=x(4)-x0
      t0=1.d0/(x1*x2*x3*(x1-x2)*(x1-x3)*(x2-x3))
      t1=x1*x2*y3
      t2=x2*x3*y1
      t3=x3*x1*y2
      t4=x1**2
      t5=x2**2
      t6=x3**2
      y1=t3*t6-t1*t5
      y3=t2*t5-t3*t4
      y2=t1*t4-t2*t6
      c1=t0*(x1*y1+x2*y2+x3*y3)
      c2=-t0*(y1+y2+y3)
      c3=t0*(t1*(x1-x2)+t2*(x2-x3)+t3*(x3-x1))
      cf(1,1)=c1
      cf(2,1)=c2
      cf(3,1)=c3
      cf(1,2)=c1+2.d0*c2*x1+3.d0*c3*t4
      cf(2,2)=c2+3.d0*c3*x1
      cf(3,2)=c3
      if (n.eq.4) then
        cf(1,3)=c1+2.d0*c2*x2+3.d0*c3*t5
        cf(2,3)=c2+3.d0*c3*x2
        cf(3,3)=c3
        cf(1,4)=c1+2.d0*c2*x3+3.d0*c3*t6
        cf(2,4)=c2+3.d0*c3*x3
        cf(3,4)=c3
        return
      end if
      do i=3,n-2
        y0=f(i)
        y1=f(i-1)-y0
        y2=f(i+1)-y0
        y3=f(i+2)-y0
        x0=x(i)
        x1=x(i-1)-x0
        x2=x(i+1)-x0
        x3=x(i+2)-x0
        t1=x1*x2*y3
        t2=x2*x3*y1
        t3=x3*x1*y2
        t0=1.d0/(x1*x2*x3*(x1-x2)*(x1-x3)*(x2-x3))
        c3=t0*(t1*(x1-x2)+t2*(x2-x3)+t3*(x3-x1))
        t4=x1**2
        t5=x2**2
        t6=x3**2
        y1=t3*t6-t1*t5
        y2=t1*t4-t2*t6
        y3=t2*t5-t3*t4
        cf(1,i)=t0*(x1*y1+x2*y2+x3*y3)
        cf(2,i)=-t0*(y1+y2+y3)
        cf(3,i)=c3
      end do
      c1=cf(1,n-2)
      c2=cf(2,n-2)
      c3=cf(3,n-2)
      cf(1,n-1)=c1+2.d0*c2*x2+3.d0*c3*t5
      cf(2,n-1)=c2+3.d0*c3*x2
      cf(3,n-1)=c3
      cf(1,n)=c1+2.d0*c2*x3+3.d0*c3*t6
      cf(2,n)=c2+3.d0*c3*x3
      cf(3,n)=c3
      return
      end subroutine


      
      
      

