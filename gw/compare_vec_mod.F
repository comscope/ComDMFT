	real*8 function compare_vec_mod(a,b,bas)
	use solid_mod
	implicit none
	real*8, intent(in) :: a(3),b(3),bas(3,3)
	integer :: i,j,k
	real*8 :: d,v(3),v1(3)
	d=100000.d0
	v1=b-a
	do k=-3,3
        do j=-3,3
	    do i=-3,3
	      v=v1+i*bas(:,1)+j*bas(:,2)+k*bas(:,3)
	      d=min(d,dot_product(v,v))
	    enddo
	  enddo
	enddo
	compare_vec_mod=d
	end