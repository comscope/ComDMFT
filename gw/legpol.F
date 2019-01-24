      subroutine legpol(p,x,lmax)
*********************************************************
*  calculate legandere polynoms up to l equal to lmax   *
*********************************************************
      implicit none
	integer :: lmax,l
      real*8 :: p(0:lmax),a,b,x
      p(0)=1.d0
      if(lmax.eq.0)return
      p(1)=x
      do l=1,lmax-1
        a=dble(2*l+1)/dble(l+1)
        b=dble(  l  )/dble(l+1)
        p(l+1)=a*x*p(l)-b*p(l-1)
	enddo
      end