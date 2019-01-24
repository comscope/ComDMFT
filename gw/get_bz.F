      subroutine get_bz
	use solid_mod
      implicit none
	integer :: i
      call cross(gbas(1,1),rbas(1,2),rbas(1,3))
      call cross(gbas(1,2),rbas(1,3),rbas(1,1))
      call cross(gbas(1,3),rbas(1,1),rbas(1,2))
      amega=rbas(1,1)*gbas(1,1)+gbas(2,1)*rbas(2,1)+gbas(3,1)*rbas(3,1)
      do i=1,3
        gbas(i,1)=gbas(i,1)/amega
        gbas(i,2)=gbas(i,2)/amega
        gbas(i,3)=gbas(i,3)/amega
      enddo
      amega=abs(amega)*par**3
      end