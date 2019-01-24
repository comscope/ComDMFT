      subroutine g_pair_print_prepare
	use atom_mod
	use manager_mod
	use parallel_mod
      use sitegen_mod
	use solid_mod
      implicit none
      integer :: i,n,ig
      real*8 :: v(3),v1(3),v2(3),diff
	allocate(g_pair_directions(3,ngroup))
	allocate(am1_g_pair_dir(ngroup))
	r_max_g_pair=par*2
	n_int_g_pair=100
	v=0.d0
	v(3)=1.d0
	g_pair_dir_number=0
      do ig=1,ngroup
        call rotate(v(1),v(2),v(3),v1(1),v1(2),v1(3),u(2,ig),1)
        n=0
        do i=1,g_pair_dir_number
          v2=g_pair_directions(:,i)-v1
          diff=dot_product(v2,v2)
          if(diff<1.d-10) then
            am1_g_pair_dir(ig)=i
            n=1
          endif
        enddo
        if(n==0) then
          g_pair_dir_number=g_pair_dir_number+1
          g_pair_directions(:,g_pair_dir_number)=v1
          am1_g_pair_dir(ig)=g_pair_dir_number
        endif
      enddo
      end