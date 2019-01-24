      subroutine ordering_all(n,a,iref)
c     Orders N elements in ascending manner
c	and gives their positions in the original N-element array
      implicit none
      integer, intent(in) :: n
      integer, intent(out) :: iref(n)
      real*8, intent(in) :: a(n)
	integer :: j,nout,k
	real*8 :: s
	integer, allocatable :: ii(:)
	allocate(ii(n))
      ii=0
	do nout=1,n
        s=1.d5
        k=1
        do j=1,n
          if(ii(j).eq.0.and.a(j)<s) then
	      s=a(j)
	      k=j
	    endif
        enddo
	  ii(k)=nout
        iref(nout)=k
	enddo
	deallocate(ii)
      end