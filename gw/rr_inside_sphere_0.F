      subroutine rr_inside_sphere_0(rrr,abas,rsph,nsph)
	use atom_mod
	use solid_mod
      implicit none
      real*8, intent(in) :: rrr(3),abas(3,3)
      integer, intent(out) :: nsph
      real*8, intent(out) :: rsph(3,8)
	integer :: i1,i2,i3
	real*8 :: v(3),v1(3),rrt(3),s,s1
	nsph=1
	rsph(:,1)=rrr
	v=rrr
	s=dot_product(v,v)
	do i3=-2,2
	  do i2=-2,2
	    do i1=-2,2
	      if(i1==0.and.i2==0.and.i3==0) cycle
	      rrt=rrr+i1*abas(:,1)+i2*abas(:,2)+i3*abas(:,3)
	      v1=rrt
	      s1=dot_product(v1,v1)
	      if(abs(s1-s)<1.d-7) then
	        nsph=nsph+1
	        rsph(:,nsph)=rrt
	      else if(s1<s) then
	        nsph=1
	        s=s1
	        rsph(:,nsph)=rrt
	      endif
	    enddo
	  enddo
	enddo  
      end