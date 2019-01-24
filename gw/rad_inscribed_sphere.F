      real*8 function rad_inscribed_sphere(a)
c ------- Calculates the radius of sphere inscribed in parallelepiped ---
	use manager_mod
	use parallel_mod
	use solid_mod
      implicit none
	real*8, intent(in) :: a(3,3)
	real*8 :: r,v(3),v_abs
	call cross(v,a(1,2),a(1,3))
	v_abs=sqrt(dot_product(v,v))
	r=abs(dot_product(v,a(:,1)))/v_abs
	call cross(v,a(1,1),a(1,3))
	v_abs=sqrt(dot_product(v,v))
	r=min(r,abs(dot_product(v,a(:,2)))/v_abs)
	call cross(v,a(1,1),a(1,2))
	v_abs=sqrt(dot_product(v,v))
	r=min(r,abs(dot_product(v,a(:,3)))/v_abs)
	rad_inscribed_sphere=r/2.d0
      end