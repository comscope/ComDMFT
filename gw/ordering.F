      subroutine ordering(n,a,m,iref)
c	 Finds among N elements M with the smallest values
c	and gives their positions in the original N-element array
      implicit none
      integer, intent(in) :: n,m
      integer, intent(out) :: iref(m)
      real*8, intent(in) :: a(n)
	integer :: j,nout,k
	real*8 :: s
	integer, allocatable :: ii(:)
	allocate(ii(n))
      ii=0
	do nout=1,m
        s=1.d5
        k=0
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
