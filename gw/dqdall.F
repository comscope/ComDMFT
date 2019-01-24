!> \brief Return \f$ \int f(x) \mathrm{d}x \f$
!>
!> Given \f$ y_i = f(x_i) \mathrm{d}x_i \f$ for \f$ 0\le i \le n \f$
!> this function returns \f$ \int f(x) \mathrm{d}x \f$.
!> The function uses different quadratures for different cases.
!> If \f$ n < 6 \f$ a specific \f$ n+1 \f$ point quadrature is chosen.
!> If \f$ n \ge 6 \f$ a specific 7 point quadrature is used for the
!> first \f$ \mathrm{mod}(n,6) \f$ points, followed by a common 7
!> point quadrature for all subsequent sets of 6 points.
!>
      real*8 function dqdall(h,y,n)
      implicit none
!>
!>    The number of data points
      integer, intent(in) :: n
!>
!>    The overall weight
      real*8, intent(in) :: h
!>
!>    The \f$ f(x_i)\mathrm{d}x_i \f$ data
      real*8, intent(in) :: y(0:n)
	integer :: n0,i
	real*8 :: fac,sum06,sum15,sum24,sum3
      dqdall=0.0d0
      if(n.eq.0) then
         dqdall=h/1.d0*(y(0))
         return
      else if(n.eq.1) then
         dqdall=h/2.d0*(y(0)+y(1))
         return
      else if(n.eq.2) then
         dqdall=h/3.d0*(y(0)+y(2)+4.d0*y(1))
         return
      else if(n.eq.3) then
         dqdall=h/8.d0*(3.d0*(y(0)+y(3))+9.d0*(y(1)+y(2)))
         return
      else if(n.eq.4) then
         dqdall=h/45.d0*(14.d0*(y(0)+y(4))+64.d0*(y(1)+y(3))+24.d0*y(2))
         return
      else if(n.eq.5) then
         dqdall=h/288.d0*(95.d0*(y(0)+y(5))+375.d0*(y(1)+y(4))+
     &                    250.d0*(y(2)+y(3)))
         return
      endif
      n0=mod(n,6)
      if(n0.eq.0) then
         dqdall=0.d0
      else if(n0.eq.1) then
         fac=h/60480.d0
         dqdall=fac*(19087.d0*y(0)+65112.d0*y(1)-46461.d0*y(2)+
     &               37504.d0*y(3)-20211.d0*y(4)+6312.d0*y(5)-
     &               863.d0*y(6))
      else if(n0.eq.2) then
         fac=h/3780.d0
         dqdall=fac*(1139.d0*y(0)+5640.d0*y(1)+33.d0*y(2)+
     &               1328.d0*y(3)-807.d0*y(4)+264.d0*y(5)-
     &               37.d0*y(6))
      else if(n0.eq.3) then
         fac=h/2240.d0
         dqdall=fac*(685.d0*y(0)+3240.d0*y(1)+1161.d0*y(2)+
     &               2176.d0*y(3)-729.d0*y(4)+216.d0*y(5)-
     &               29.d0*y(6))
      else if(n0.eq.4) then
         fac=h/945.d0
         dqdall=fac*(286.d0*y(0)+1392.d0*y(1)+384.d0*y(2)+
     &               1504.d0*y(3)+174.d0*y(4)+48.d0*y(5)-
     &               8.d0*y(6))
      else if(n0.eq.5) then
         fac=h/12096.d0
         dqdall=fac*(3715.d0*y(0)+17400.d0*y(1)+6375.d0*y(2)+
     &               16000.d0*y(3)+11625.d0*y(4)+5640.d0*y(5)-
     &               275.d0*y(6))
      endif
      fac=h/140.d0
      sum06=0.d0
      sum15=0.d0
      sum24=0.d0
      sum3=0.d0
      do i=n0,n-6,6
         sum06=sum06+y(i)+y(i+6)
         sum15=sum15+y(i+1)+y(i+5)
         sum24=sum24+y(i+2)+y(i+4)
         sum3=sum3+y(i+3)
      enddo
      dqdall=dqdall+fac*(41.d0*sum06+216.d0*sum15+27.d0*sum24+
     &                   272.d0*sum3)
      end
      
      subroutine dqd(h,y,z,ndim)
*****************************************************************
*       the subroutine calculates ndim integrals                *
*       between x1 and xn of function z ( n.le.ndim             *
*        xn=x1+(n-1)*h ) for the special case:                  *
*       ndim=2*nd+1 y(1)=0 . z(2*i+1) is calculated             *
*       by the simpson's method and z(2*i+2)=z(2*i+1)+		*
*       delt(2*i+2) . to increase a accuracy it is              *
*       used delt(2*i+2)=delt4(2*i+2)*(z(2*i+3)-z(2*i+1))/      *
*       (delt4(2*i+2)+delt4(2*i+3)) where delt4(2*i+2)          *
*       is calculated by using a cubic interpolation            *
*       between y(2*i),y(2*i+1),y(2*i+2),y(2*i+3).              *
* Revised by Yu.A.Uspensky (FIAN)                               *
*****************************************************************
c Yurii Alekseevich Uspenskii, P.N. Lebedev Institute, Moscow.
c FIAN == Physical Institute of the Academy of Sciences
c      == P.N. Lebedev Institute
	use manager_mod
	use parallel_mod
      implicit none
	integer, intent(in) :: ndim
      real*8, intent(in) :: h,y(ndim)
      real*8, intent(out) :: z(ndim)
	integer :: i
	real*8 :: c3,c24,y2,y3,y4,y5,sum2,y1,delt3,delt4,delt5,delt
c
      if((ndim/2+ndim/2).lt.ndim) go to 1
        if(maswrk) write(iun,*)'the subroutine qd works only if
     &    ndim=2*nd+1'
        call ending
c
1     c3=h/3.d0
      c24=h/24.d0
      z(1)=0.d0
      if(ndim.eq.1) go to 101
      y2=y(1)
      y3=y(2)
      y4=y(3)
      z(2)=c24*(10*y2+16*y3-2*y4)
c       z(2)=(y2+y3)/2.d0*h
      sum2=y2+4*y3+y4
      z(3)=c3*sum2
      if(ndim.eq.3) go to 101
      y5=y(4)
c
c       the main loop of a integration
c
      do i=4,ndim-1,2
        y1=y3
        y2=y4
        y3=y5
        y4=y(i+1)
        delt3=y2+4*y3+y4
        sum2=sum2+delt3
        z(i+1)=c3*sum2
        delt4=-y1+13*y2+13*y3-y4
        if(i.ne.ndim-1) go to 3
        delt5=delt3*8
        go to 4
3       y5=y(i+2)
        delt5=-y1+12*y2+26*y3+12*y4-y5
4       continue
        if(abs(delt5).gt.1.e-20) then
          delt=8*delt4*delt3/delt5
          z(i)=z(i-1)+c24*delt
        else
          z(i)=z(i-1)+c24*delt4
        end if
	enddo
101   end
