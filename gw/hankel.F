      subroutine hankel(r,e,s,han,lmax)
************************************************************
*  calculate hankel functions h defined as follows:        *
*                     (ks)^l+1                             *
*            h (kr)=- -------- h (kr)                      *
*             l       (2l-1)!!  l                          *
*  and their derevatives with respect to r. here k=sqrt(e) *
*  and h(x) is a spherical hankel function which can be    *
*  expressed through the neuman and bessel functions       *
*                h(x)=n(x)-ij(x)                           *
************************************************************
	implicit none
	integer, intent(in) :: lmax
	real*8, intent(in) :: r,e,s
	real*8, intent(out) :: han(3,0:lmax)
	integer :: l
	real*8 :: x,t,hant(3,0:1)
	x=r*sqrt(-e)
	t=s*sqrt(-e)
	hant(1,0)=t/x*exp(-x)
	hant(1,1)=hant(1,0)*t*(1.d0+1.d0/x)
	hant(2,0)=-1.d0/s*hant(1,1)
	hant(2,1)=-2.d0/r*hant(1,1)-t*t/s*hant(1,0)
	hant(3,0)=-e*hant(1,0)-2.d0/r*hant(2,0)
	hant(3,1)=(2.d0/r/r-e)*hant(1,1)-2.d0/r*hant(2,1)
	do l=0,min(lmax,1)
	  han(1,l)=hant(1,l)
	  han(2,l)=hant(2,l)
	  han(3,l)=hant(3,l)
	enddo
	do l=1,lmax-1
	  han(1,l+1)=t/x*han(1,l)+t*t/dble((2*l-1)*(2*l+1))*han(1,l-1)
	  han(2,l+1)=-dble(l+2)/r*han(1,l+1)-t*t/s/dble(2*l+1)*han(1,l)
	  han(3,l+1)=(dble((l+1)*(l+2))/r/r-e)*han(1,l+1)-
     &             2.d0/r*han(2,l+1)
      enddo
	end