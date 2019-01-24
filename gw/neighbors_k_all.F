      subroutine neighbors_k_all(k,qsph,nsph)
	use atom_mod
	use solid_mod
      implicit none
      integer, intent(in) :: k
      integer, intent(out) :: nsph
      real*8, intent(out) :: qsph(3,12)
	integer :: i1,i2,i3,kc
	real*8 :: v(3),v1(3),s,s1
	nsph=0
	v=pnt(:,k)
	s=1.d5
	do kc=1,nqdiv
	  do i3=-2,2
	    do i2=-2,2
	      do i1=-2,2
	        v1=pnt(:,kc)+i1*gbas(:,1)+i2*gbas(:,2)+i3*gbas(:,3)
	        v1=v1-v
	        s1=dot_product(v1,v1)
	        if(abs(s1-s)<1.d-7) then
	          nsph=nsph+1
	          qsph(:,nsph)=v1+v
	        else if(s1<s) then
	          nsph=1
	          s=s1
	          qsph(:,nsph)=v1+v
	        endif 
	      enddo
	    enddo
	  enddo 
	enddo
      end