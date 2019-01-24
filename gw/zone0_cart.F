C-----------------------------------------------------------------------
C
C       Map the translation vectors into the Brillouin zone
C
C-----------------------------------------------------------------------
C
C       Take the translation vector associated with a symmetry 
C       operation and add different multiples of the lattice vectors.
C       Check which translation is the shortest overall and select
C       that one.
C
C       Comments:
C       1. In practice q and q1 are the same variables when this 
C          routine is called in subroutine group. This is potentially
C          dangerous as the Fortran standard assumes that all arguments
C          are non-overlapping pieces of memory. Hence compiler 
C          optimizations may break code if this assumption is false.
C          To avoid problems the original data is kept in v2 so that
C          input data may be overwritten without causing problems.
C       2. Different translation vectors constructed may be equally
C          long, for example a vector 0.5 is equally long as the 
C          vector 0.5-1=-0.5. The simple comparison d0<d may then pick
C          a different result depending on small deviations from the
C          correct value. For example 0.5+delta is larger than 
C          -0.5+delta. To avoid this we look for vectors that are
C          significantly shorter by comparing d0+1.0d-6<d.
C       3. In practice the code seems to work better if the selected
C          translation vector is aligned with the lattice vectors. 
C          For this reason we want to try those vectors first. Therefore
C          we run the i, j, and k loops in reverse order, i.e. from 3
C          to -3.
C
	subroutine zone0_cart(q,bas,q1)
	use solid_mod
	implicit none
	real*8, intent(in) :: q(3),bas(3,3)
	real*8, intent(out) :: q1(3)
	integer :: i,j,k
	real*8 :: d0,d,v(3),v1(3),v2(3)
	v1=q
	v2=q
        d=dot_product(v1,v1)
	do k=3,-3,-1
          do j=3,-3,-1
	    do i=3,-3,-1
	      v=v2+i*bas(:,1)+j*bas(:,2)+k*bas(:,3)
	      d0=dot_product(v,v)
	      if(d0+1.0d-6<d) then
	        d=d0
	        v1=v
	      endif
	    enddo
	  enddo
	enddo
	q1=v1
	end
C
C-----------------------------------------------------------------------
