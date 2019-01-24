      subroutine gtbvec(k,b,shift,v)
      implicit none
	integer, intent(in) :: k(3)
      real*8, intent(in) :: b(3,3),shift(3)
      real*8,intent(out) :: v(3)
      v(1) = shift(1) + k(1)*b(1,1) + k(2)*b(1,2) + k(3)*b(1,3)
      v(2) = shift(2) + k(1)*b(2,1) + k(2)*b(2,2) + k(3)*b(2,3)
      v(3) = shift(3) + k(1)*b(3,1) + k(2)*b(3,2) + k(3)*b(3,3)
      end