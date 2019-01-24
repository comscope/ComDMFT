      subroutine legendre(n,x,p)
c- calculates legendre polynomical using a recursion relation
c ----------------------------------------------------------------------
ci inputs
ci   n,x
co outputs
co   p(0:n;x)
cr remarks
cr   recursion relation is
cr   p_n = ((2*n-1)*x*p_(n-1) - (n-1)*p_(n-2))/n
c ----------------------------------------------------------------------
      implicit none
c passed parameters
      integer, intent(in) :: n
      real*8, intent(in) :: x
      real*8, intent(out) :: p(0:n)
c local parameters
      integer :: l,n1
      p(0)=1.d0
      if(n==0) return
      p(1)=x
      if(n==1) return
      do l=2,n
        n1=l-1
        p(l)=((n1+l)*x*p(l-1)-n1*p(l-2))/l
	enddo
      end