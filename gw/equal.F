      subroutine equal(key,a,eps,n)
c     Makes the elements of array A exactly equal when they differ
c                         less than by EPS      
c     key =  1 - for Real array
c     key =  2 - for Complex array
	use atom_mod
	use etot_mod
	use manager_mod
	use solid_mod
	use units_mod
      implicit none
      integer, intent(in) :: key,n
      real*8, intent(in) :: eps
      real*8, intent(inout) :: a(n,key)
	integer :: i,j
	real*8 :: d
	integer, allocatable :: ind(:)
	allocate(ind(n))
	ind=0
	do i=1,n
	  if(ind(i)==1) cycle
	  do j=i+1,n
	    d=abs(a(j,1)-a(i,1))
	    if(key==2) d=d+abs(a(j,2)-a(i,2))
	    if(d<eps) then
	      a(j,1:key)=a(i,1:key)
	      ind(j)=1
	    endif
	  enddo
	  ind(i)=1
	enddo
	deallocate(ind)
      end