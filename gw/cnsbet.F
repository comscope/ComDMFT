        subroutine cnsbet(ils,bett)
***************************************************************
*                                                             *
*   program for calculating the norm. coefficients bett       *
*                                                             *
*               / (2l+1)*(l-n)! \ (1/2)                       *
*    bett(l,m)=(  -------------  )      ; n=abs(m)            *
*               \  4*pi *(l+n)! /                             *
*                                                             *
*  written by s.yu.savrasov (p.n.lebedev physical institute)  *
***************************************************************
	use units_mod
      implicit none
	integer, intent(in) :: ils
      real*8, intent(out) :: bett((ils+1)*(ils+1))
	integer :: l,m,lm
	real*8 :: f,fact12
      do l=0,ils
        do m=-l,l
          lm=l*(l+1)+m+1
          f=fact12(l-abs(m),l+abs(m))
          if(m.ne.0) f=2.d0*f
          bett(lm)=sqrt((2.d0*l+1.d0)/4.d0/pi*f)
	  enddo
	enddo
      end