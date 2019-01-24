      subroutine rotate_opt(v,v1,u,key)
	use units_mod
      implicit none
	integer,intent(in) :: key
      real*8, intent(in) :: u(3,3),v(3)
      real*8, intent(out) :: v1(3)
      if(key.eq.1) v1=matmul(v,u) !!! inversional ig-operation
      if(key.eq.2) v1=matmul(u,v) !!!      direct ig-operation
      end