      subroutine group_table
	use manager_mod
	use parallel_mod
	use solid_mod
      implicit none
	integer :: ig,igg,iggg
	real*8 :: d(3,3),a(3),dif,compare_vec_mod
	allocate(u_prod(ngroup,ngroup))
	allocate(u_inv(ngroup))
c ------ Inverse elements -----------------------------------------	
	do ig=1,ngroup
	  d=transpose(u_opt(:,:,ig))
	  a=-matmul(d,shift(:,ig))
	  do igg=1,ngroup
	    dif=maxval(abs(d-u_opt(:,:,igg)))
     &	     +compare_vec_mod(a,shift(1,igg),rbas)
	    if(dif<1.d-8) u_inv(ig)=igg
	  enddo
	enddo
c ------ Multiplication table -----------------------------------------	
	do igg=1,ngroup
	  do ig=1,ngroup
	    d=matmul(u_opt(:,:,ig),u_opt(:,:,igg))
	    a=matmul(u_opt(:,:,ig),shift(:,igg))+shift(:,ig)
	    do iggg=1,ngroup
	      dif=maxval(abs(d-u_opt(:,:,iggg)))
     &	       +compare_vec_mod(a,shift(1,iggg),rbas)
	      if(dif<1.d-8) u_prod(ig,igg)=iggg
	    enddo
	  enddo
	enddo
      end