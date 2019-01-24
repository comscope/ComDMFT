	subroutine set_pairs
	use solid_mod
	integer :: i,j,k
	real*8 :: x(3),d(3)
	allocate(ind_dif_tt(natom,natom))
	allocate(dif_tt(3,natom*natom))
	num_dif_tt=0
	do j=1,natom
	  do i=1,natom
	    x=tau(:,j)-tau(:,i)
	    do k=1,num_dif_tt
	      d=dif_tt(:,k)-x
	      if(abs(d(1))<1.d-9.and.abs(d(2))<1.d-9
     &		                  .and.abs(d(3))<1.d-9) then
	        ind_dif_tt(i,j)=k
              goto 1
	      endif
	    enddo
	    num_dif_tt=num_dif_tt+1
	    dif_tt(:,num_dif_tt)=x
	    ind_dif_tt(i,j)=num_dif_tt
1         continue
        enddo
	enddo
	end