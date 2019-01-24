      integer function rbs_number(v)
c     Gives the position of vector V in the list of RBS-vectors -----      
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	real*8, intent(in) :: v(3)
	integer :: i1,i2,i3
	real*8 :: det,det3
      det=det3(rbas(1,1),rbas(1,2),rbas(1,3))
      i1=idnint(det3(v,rbas(1,2),rbas(1,3))/det)
      i2=idnint(det3(rbas(1,1),v,rbas(1,3))/det)
      i3=idnint(det3(rbas(1,1),rbas(1,2),v)/det)
	rbs_number=ind_rsp(i1,i2,i3)
      end