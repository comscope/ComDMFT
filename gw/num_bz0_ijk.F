	integer function num_bz0_ijk(i,j,k)
c     Gives number of equivalent k-point in zero's BZ,
c                where (i,j,k) - indexes in first BZ	
	use solid_mod
	implicit none
	integer, intent(in) :: i,j,k
	integer :: num_bz1_ijk
	num_bz0_ijk=index_k1(num_bz1_ijk(i,j,k))
	end
	integer function num_bz1_ijk(i,j,k)
c     Gives number of equivalent k-point in first BZ,
c                where (i,j,k) - indexes in first BZ	
	use solid_mod
	implicit none
	integer, intent(in) :: i,j,k
	integer :: i0,j0,k0
	i0=mod(i,ndiv(1))
	if(i0.lt.0) i0=i0+ndiv(1)
	j0=mod(j,ndiv(2))
	if(j0.lt.0) j0=j0+ndiv(2)
	k0=mod(k,ndiv(3))
	if(k0.lt.0) k0=k0+ndiv(3)
	num_bz1_ijk=ndiv(1)*ndiv(2)*k0+ndiv(1)*j0+i0+1
	end