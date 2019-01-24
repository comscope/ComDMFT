	subroutine add_to_output_w(ind_k,k0,n_pbt,w_c_nu,key)
c	key = 0 Polarizability
c	key = 1 W_full - V
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer, intent(in) :: ind_k,k0,n_pbt,key
  	real*8, intent(in) :: w_c_nu(n_pbt,n_pbt,ndim3_nu)
	integer :: ind_nu,i_nu
	complex*16 :: wint
	complex*16, allocatable :: work(:,:)
	allocate(work(n_pbt,n_pbt))
	do ind_nu=1,ndim3_nu
	  i_nu=me3_tau*ndim3_nu+ind_nu-1
	  call unpack_hermit(work,w_c_nu(1,1,ind_nu),n_pbt,n_pbt,n_pbt,
     &                     (0.d0,0.d0),(1.d0,0.d0))
	  if(key==0) then
          call pw_from_pb(ind_k,n_pbt,n_pbt,1,work,wint,3)
	    p_output(i_nu,k0)=real(wint)
	  else if(key==1) then
          call pw_from_pb(ind_k,n_pbt,n_pbt,1,work,wint,1)
	    w_output(i_nu,k0)=real(wint)
        endif
	enddo
	deallocate(work)
      end