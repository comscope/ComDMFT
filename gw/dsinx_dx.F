      real*8 function dsinx_dx(z)
c We use series expansion, so z should not be too large --------------      
      implicit none
      real*8, intent(in) :: z
      integer :: k,k2
      real*8 :: f,s,z2
      s=-z/6.d0
      f=s
      z2=z*z
      do k=2,10000
        k2=k+k
        f=-f*z2/2.d0/(k2+1.d0)
        s=s+f
        if(abs(f)<1.d-30) exit
      enddo
      dsinx_dx=s
      end