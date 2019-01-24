      subroutine s_overlap_m_b(s,n,n1,iref,iref1)
c     Calculates interstitial overlap matrix for plane waves ------
c     First plane wave (left) belongs to the basis with its IREF
c     Second plane wave (right) belongs to the basis with its IREF1
	use atom_mod
	use etot_mod
	use manager_mod
	use solid_mod
	use units_mod
      implicit none
      integer, intent(in) :: n,n1,iref(n),iref1(n1)
      complex*16, intent(out) :: s(n,n1)
	integer :: j,j0,i,i0,ind,gbs_number
	real*8 :: v(3)
	do j=1,n1
	  j0=iref1(j)
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
      end