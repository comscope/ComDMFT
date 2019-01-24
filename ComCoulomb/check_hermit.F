      double precision function check_hermit(a,n,n0)
      implicit none
      integer :: n,i,j,n0
      double precision :: d
      complex*16 :: a(n0,n0)
      d=0.d0
      do j=1,n
        do i=1,j
          d=max(d,abs(a(i,j)-conjg(a(j,i))))
          a(i,j)=0.5d0*(a(i,j)+conjg(a(j,i)))
          a(j,i)=conjg(a(i,j))
        enddo
      enddo
      check_hermit=d
      end
      
