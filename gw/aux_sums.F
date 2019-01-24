      double precision function csum0(x)
c     sums (-1)^n x^(2n)/(2n+1)! for n=0 to inf
      implicit none
      real*8, intent(in) :: x
	integer :: n
	real*8 :: x2,s,f
	x2=x*x
	s=1.d0
	f=s
	n=0
1     n=n+1
      f=-f*x2/2/n/(n+n+1)
      s=s+f
      if(abs(f)>1.d-30) goto 1
      csum0=s
      end
      
      double precision function ssum1(x)
c     sums (-1)^n x^(2n)/(2n+1)!/(2n+3) for n=0 to inf
      implicit none
      real*8, intent(in) :: x
	integer :: n
	real*8 :: x2,s,f
	x2=x*x
	f=1.d0
	s=f/3.d0
	n=0
1     n=n+1
      f=-f*x2/2/n/(n+n+1)
      s=s+f/(n+n+3)
      if(abs(f)>1.d-30) goto 1
      ssum1=s
      end
      
      double precision function csum2(x)
c     sums (-1)^n x^(2n-2)*(2n-3)*n/(2n+1)! for n=2 to inf
      implicit none
      real*8, intent(in) :: x
	integer :: n
	real*8 :: x2,s,f
	x2=x*x
	f=x2/120.d0
	s=f+f
	n=2
1     n=n+1
      f=-f*x2/2/n/(n+n+1)
      s=s+f*n*(n+n-3)
      if(abs(f)>1.d-30) goto 1
      csum2=s
      end
      
      double precision function ssum3(x)
c     sums (-1)^n x^(2n-3)*(2n^2-9n+7)/(4n^2-1)! for n=2 to inf
      implicit none
      real*8, intent(in) :: x
	integer :: n
	real*8 :: x2,s,f,a
	x2=x*x
	f=x/2.d0
	s=-f/5.d0
	n=2
1     n=n+1
      a=n+n-2
      f=-f*x2/a/(a-1.d0)
      s=s+f*real(2*n*n-9*n+7)/(4*n*n-1)
      if(abs(f)>1.d-30) goto 1
      ssum3=s
      end