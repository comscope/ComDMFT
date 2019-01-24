      subroutine int1_xn_exp_ax(n,a,b,h,e,de,t)    
      implicit none
      integer, intent(in) :: n
      real*8, intent(in) :: a,b,h,de,t
      real*8, intent(out) :: e(0:n)
      integer :: i
      real*8 :: eca,ecb,green0,c
      c=de*h
      eca=green0(de,t-h*a)
      ecb=green0(de,t-h*b)
      e(0)=h*(ecb-eca)/c
      do i=1,n
        e(i)=(h*(b**i*ecb-a**i*eca)-i*e(i-1))/c
      enddo
      end