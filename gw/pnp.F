      real*8 function pnp(n,x)
c- calculates derivative of legendre polynomical from recursion relation
c ----------------------------------------------------------------------
ci inputs
ci   n,x
co outputs
co   pnp: p'_n(x)
cr remarks
cr   recursion relations for p and p' are
cr   p_n (x) = ((2*n-1)*x*p_(n-1) - (n-1)*p_(n-2))/n
cr   p'_n(x) = n/(1-x^2) (-x*p_n + p_(n-1) )
c ----------------------------------------------------------------------
      implicit none
c passed parameters
      integer n
      real*8 :: x
c local parameters
      real*8 :: jpjm1,cj,pjp1,ln
      integer j

c jpjm1 is j*p_(j-1);  cj is 2*j - 1;  pjp1 is p_(j+1)
      jpjm1 = 0
      ln = 1
      cj = 1
      do j = 1, n
        pjp1 = (cj*x*ln - jpjm1)/j
        jpjm1 = j*ln
        cj = cj + 2
        ln = pjp1
	enddo
      pnp = 1/(1-x**2)*(-n*x*ln + jpjm1)
      end