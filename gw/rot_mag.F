      subroutine rot_mag(ig,b_old,b_new)
	use solid_mod
      implicit none
	integer :: ig,i,j
	real*8 :: b_new(3),b_old(3)
      do i=1,3
	  b_new(i)=0.d0
        do j=1,3
          b_new(i)=b_new(i)+agm(i,j,ig)*b_old(j)
        enddo  !!! over j
	  b_new(i)=dfloat(ifia(ig))*b_new(i)
      enddo  !!! over i
      end