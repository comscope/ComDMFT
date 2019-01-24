      real*8 function plegn(n,x)
c- calculates legendre polynomical using a recursion relation
c ----------------------------------------------------------------------
ci inputs
ci   n,x
co outputs
co   plegn: p_n(x)
cr remarks
cr   recursion relation is
cr   p_n = ((2*n-1)*x*p_(n-1) - (n-1)*p_(n-2))/n
c ----------------------------------------------------------------------
      implicit none
c passed parameters
      integer :: n
      real*8 :: x
c local parameters
      real*8 :: jpjm1,cj,pjp1
      integer :: j

c jpjm1 is j*p_(j-1);  cj is 2*j - 1;  pjp1 is p_(j+1)
      jpjm1 = 0
      plegn = 1
      cj = 1
      do j = 1, n
        pjp1 = (cj*x*plegn - jpjm1)/j
        jpjm1 = j*plegn
        cj = cj + 2
        plegn = pjp1
	enddo
      end