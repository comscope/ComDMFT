      subroutine s_overlap(key,s,n0,n,iref)
c     Calculates interstitial overlap matrix for plane waves ------
c     IREF - reference of the plane waves to GBS array ------------
c     key = 0  Just overlap matrix
c     key = 1  Its inverse
	use atom_mod
	use etot_mod
	use manager_mod
	use solid_mod
	use units_mod
      implicit none
      integer, intent(in) :: key,n,n0,iref(n0)
      complex*16, intent(out) :: s(n0,n0)
	integer :: j,j0,i,i0,ind,info,gbs_number
	real*8 :: v(3)
	integer, allocatable :: ipiv(:)
	complex*16, allocatable :: work(:)
	do j=1,n
	  j0=iref(j)
	  do i=1,n
	    i0=iref(i)
	    v=gbs(:,j0)-gbs(:,i0)
          ind=gbs_number(v)
	    if(complex_ro) then
		    s(i,j)=dcmplx(sovr(ind),sovi(ind))
	    else
		    s(i,j)=dcmplx(sovr(ind),0.d0)
	    endif
	  enddo
	enddo
	s=s*amega
	if(key==1) then
	  allocate(ipiv(n),work(3*n))
	  call zgetrf(n,n,s,n0,ipiv,info)
	  call zgetri(n,s,n0,ipiv,work,3*n,info)
	  deallocate(ipiv,work)
      endif
      end