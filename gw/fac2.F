      real*8 function fac2(l)
	integer :: l,i,ll
      if(l.eq.-1) then
        fac2=1
        return
      endif
      if(l.eq.1) then
        fac2=1
        return
      endif
      fac2=1
      ll=(l-1)/2
      do i=1,ll
        fac2=fac2*(2*i+1)
	enddo
      end
      
      subroutine facc(f,n)
      implicit none
	integer, intent(in) :: n
      real*8, intent(out) :: f(0:n)
	integer :: i
      f(0)=1.d0
      do i=1,n
        f(i)=f(i-1)*i
      enddo
      end
      
      real*8 function fact12(n,m)   !!!   n!/m!
c***************************************************
      implicit none
	integer, intent(in) :: n,m
	integer :: i
      fact12=1.0d0
      if(n.lt.m)then
        fact12=1.d0/dble(n+1)
        do i=n+2,m
          fact12=fact12/dble(i)
        enddo
!     elseif(n.eq.m)then
!       fact12=1.d0
      elseif(n.gt.m)then
        fact12=dble(m+1)
        do i=m+2,n
          fact12=fact12*dble(i)
        enddo
      endif
      end
