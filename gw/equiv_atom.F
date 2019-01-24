	subroutine equiv_atom(t,iatom0,r0)
	use solid_mod
	implicit none
	real*8, intent(in) :: t(3)
	integer, intent(out) :: iatom0
	real*8, intent(out) :: r0(3)
	integer :: i,j,k,iatom
	real*8 :: d0,v(3),v1(3)
	do k=-3,3
        do j=-3,3
	    do i=-3,3
	      v=t+i*rbas(:,1)+j*rbas(:,2)+k*rbas(:,3)
	      do iatom=1,natom
	        v1=v-tau(:,iatom)
	        d0=dot_product(v1,v1)
	        if(d0<1.d-6) then
	          iatom0=iatom
	          r0=tau(:,iatom0)-t
	        endif
	      enddo
	    enddo
	  enddo
	enddo
	end