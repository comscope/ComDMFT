      real*8 function sin_integral(z)
c We use series expansion, so z should not be too large --------------      
      implicit none
      real*8, intent(in) :: z
      integer :: k,k2
      real*8 :: f,s
      s=z
      f=z
      do k=1,10000
        k2=k+k
        f=-f*(z/(k2+1.d0))**2*(k2-1.d0)/k2
        s=s+f
        if(abs(f)<1.d-30) exit
      enddo
      sin_integral=s
      end