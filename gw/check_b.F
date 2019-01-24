      subroutine check_b
	use manager_mod
	use parallel_mod
	use solid_mod
      implicit none
	integer :: ig,i
	real*8 :: b_new(3),dy
      do ig=1,ngroup
        ag(1,1,ig)=u(10,ig)
        ag(1,2,ig)=u(4,ig)
        ag(1,3,ig)=u(7,ig)
        ag(2,1,ig)=u(8,ig)
        ag(2,2,ig)=u(2,ig)
        ag(2,3,ig)=u(5,ig)
        ag(3,1,ig)=u(9,ig)
        ag(3,2,ig)=u(3,ig)
        ag(3,3,ig)=u(6,ig)
        agm(1,1,ig)=ag(2,2,ig)*ag(3,3,ig)-ag(3,2,ig)*ag(2,3,ig)
        agm(1,2,ig)=ag(2,3,ig)*ag(3,1,ig)-ag(3,3,ig)*ag(2,1,ig)
        agm(1,3,ig)=ag(2,1,ig)*ag(3,2,ig)-ag(3,1,ig)*ag(2,2,ig)
        agm(2,1,ig)=ag(3,2,ig)*ag(1,3,ig)-ag(1,2,ig)*ag(3,3,ig)
        agm(2,2,ig)=ag(3,3,ig)*ag(1,1,ig)-ag(1,3,ig)*ag(3,1,ig)
        agm(2,3,ig)=ag(3,1,ig)*ag(1,2,ig)-ag(1,1,ig)*ag(3,2,ig)
        agm(3,1,ig)=ag(1,2,ig)*ag(2,3,ig)-ag(2,2,ig)*ag(1,3,ig)
        agm(3,2,ig)=ag(1,3,ig)*ag(2,1,ig)-ag(2,3,ig)*ag(1,1,ig)
        agm(3,3,ig)=ag(1,1,ig)*ag(2,2,ig)-ag(2,1,ig)*ag(1,2,ig)
	  call rot_mag(ig,b_ext,b_new)
	  dy=0.d0
        do i=1,3
	    dy=dy+abs(b_new(i)-b_ext(i))
        enddo  !!! over i
	  if(dy.gt.1.d-6) then
          if(maswrk) then
            write(iun,*)'B_EXT is not invariant for ig =',ig
          endif
          call ending
	  endif
      enddo  !!! over ig
      end